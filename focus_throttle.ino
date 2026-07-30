// focus_throttle.ino - thermloop, Arduino UNO Rev3
//
// Closes the thermloop platform around a live BCI workload signal.
//
//   OpenBCI GUI (host)  --Serial 57600-->  UNO
//        Focus Widget                       |
//        [0,1] per line                     |
//                                           v
//                            focus -> requested power P_want
//                                           |
//                            ThermalThrottle::solve()  -> alpha
//                                           |
//                    alpha * P_want --> PWM D6 --> power resistor
//                                           |
//                            NTC on A0 --> measured T --> syncTemperature()
//                                           |
//                            PIController --> 25 kHz PWM D9 --> fan
//
// The host does the EEG DSP; the UNO never sees a sample of EEG. Focus
// arrives as a workload DEMAND, which is exactly the planned_P input the
// predictive throttle was built to take.
//
// Wiring:
//   A0  NTC divider (10k NTC to GND, 10k fixed to 5V)
//   D6  MOSFET gate driving the power resistor
//   D9  fan PWM  (OC1A -- reconfigured to 25 kHz, see setupFanPWM)
//
// Serial protocol matches the OpenBCI reference sketch: 57600 baud,
// newline terminated, 32 byte receive buffer. Send "CAL" to run a step
// response for identify.cpp.

#include "throttle.h"
#include "controller.h"

// ---------------------------------------------------------------------
// Pins and protocol
// ---------------------------------------------------------------------
const byte PIN_NTC  = A0;
const byte PIN_HEAT = 6;
const byte PIN_FAN  = 9;

const long BAUD     = 57600;
const byte numChars = 32;
char receivedChars[numChars];
boolean newData = false;

// ---------------------------------------------------------------------
// Control timing
// ---------------------------------------------------------------------
const float DT       = 0.5;    // seconds per control tick
const int   HORIZON  = 24;     // coarse steps looked ahead
const int   STRIDE   = 5;      // 5 * 0.5 s = 2.5 s per coarse step
                               // -> 60 s of lookahead for 24 iterations
const float T_LIMIT  = 45.0;   // C, hard ceiling on the rig
const float T_TARGET = 38.0;   // C, fan setpoint (below the ceiling)
const float P_MAX    = 2.0;    // W, resistor at full duty

// ---------------------------------------------------------------------
// Plant.
//
// These came from a step response fitted by identify.cpp, NOT from a
// plausible-looking guess. Re-run the identification whenever the rig
// changes -- new sink, new bonding, new enclosure -- because every
// prediction in this sketch rests on them.
//
//   1. Send "CAL" and capture the CSV that follows.
//   2. ./identify cal.csv 2.0
//   3. Paste the two lines it prints here.
//
// The values below are the self-test defaults and WILL be wrong for your
// hardware; treat them as a shape, not a measurement.
// ---------------------------------------------------------------------
const float PLANT_R[2] = {4.1474, 8.8781};    // C/W
const float PLANT_C[2] = {1.5880, 14.6102};   // J/C
const float T_AMB      = 24.0;                // C, room reference

// Sync time constant: several times the slowest plant tau (about 130 s),
// so the offset estimator tracks ambient drift without chasing NTC noise.
const float SYNC_TAU = 600.0;

ThermalThrottle throttle;
PIController fan(0.15, 0.02, DT, 0.0, 1.0);   // clamped: see controller.h

// ---------------------------------------------------------------------
// Workload plan
//
// The Focus Widget emits a binary concentration flag; smoothing gives a
// continuous demand level. What the plan does with that level across the
// horizon is a real modelling choice:
//
//   envelope = 0  pure persistence. Assumes the current level holds.
//                 Cheap and responsive, but blind to a focus spike that
//                 arrives mid-horizon.
//   envelope = 1  rising envelope. Demand ramps from the current level to
//                 full burst by the end of the horizon, so the throttle
//                 budgets for the worst case the user could plausibly
//                 reach. Costs some headroom; never caught out.
//
// Default is 1. A thermal ceiling is a constraint, not a target, and the
// re-solve every tick recovers any headroom the envelope gave away. The
// asymmetry decides it: over-throttling costs a little responsiveness,
// under-throttling costs an overshoot past T_LIMIT.
// ---------------------------------------------------------------------
struct FocusPlan {
    float idle;        // W, telemetry floor
    float burstNow;    // W, burst at the current smoothed focus
    float burstMax;    // W, burst at full focus
    int   period;      // coarse steps per burst cycle
    int   duty;        // coarse steps of burst within a cycle
    int   horizon;     // coarse steps in the plan
    float envelope;    // 0 = persistence, 1 = ramp to burstMax
};

FocusPlan plan = {0.10, 0.30, 1.60, 4, 2, HORIZON, 1.0};

float planAt(int k, const void* ctx) {
    const FocusPlan* p = (const FocusPlan*)ctx;
    if ((k % p->period) >= p->duty) return p->idle;
    float ramp = p->envelope * (float)k / (float)p->horizon;
    if (ramp > 1.0) ramp = 1.0;
    return p->burstNow + (p->burstMax - p->burstNow) * ramp;
}

float focusLevel = 0.0;   // EMA of the binary focus flag, 0..1

// ---------------------------------------------------------------------
// Fan PWM at 25 kHz.
//
// A 4-wire fan expects 25 kHz on its control line; the Arduino default of
// roughly 490 Hz sits inside the audible band and many fans simply ignore
// it. Timer1 in Fast PWM with ICR1 as TOP gives an exact rate:
//   f = 16 MHz / (prescale * (1 + TOP)) = 16e6 / (1 * 640) = 25 kHz
//
// This takes Timer1 away from analogWrite() on D9 and D10 -- write OCR1A
// directly instead. millis() lives on Timer0 and is unaffected.
// ---------------------------------------------------------------------
const unsigned int FAN_TOP = 639;

void setupFanPWM() {
    pinMode(PIN_FAN, OUTPUT);
    TCCR1A = _BV(COM1A1) | _BV(WGM11);              // non-inverting, WGM11
    TCCR1B = _BV(WGM13) | _BV(WGM12) | _BV(CS10);   // ICR1 top, prescale 1
    ICR1 = FAN_TOP;
    OCR1A = 0;
}

void setFan(float duty) {
    if (duty < 0.0) duty = 0.0;
    if (duty > 1.0) duty = 1.0;
    OCR1A = (unsigned int)(duty * (float)FAN_TOP);
}

void setHeater(float duty) {
    if (duty < 0.0) duty = 0.0;
    if (duty > 1.0) duty = 1.0;
    analogWrite(PIN_HEAT, (int)(duty * 255.0));
}

// ---------------------------------------------------------------------
// NTC read. Beta model on a 10k/10k divider.
// Get beta from your datasheet, or from identify.cpp's two-point solver.
// ---------------------------------------------------------------------
const float NTC_R25  = 10000.0;
const float NTC_BETA = 3950.0;
const float R_FIXED  = 10000.0;

float readTemperature() {
    int raw = analogRead(PIN_NTC);
    if (raw <= 0) raw = 1;              // guard an open or shorted divider
    if (raw >= 1023) raw = 1022;
    float r = R_FIXED * (float)raw / (1023.0 - (float)raw);
    float invT = 1.0 / 298.15 + log(r / NTC_R25) / NTC_BETA;
    return 1.0 / invT - 273.15;
}

// ---------------------------------------------------------------------
// Serial receive, after the OpenBCI reference sketch
// ---------------------------------------------------------------------
void recvWithEndMarker() {
    static byte ndx = 0;
    char endMarker = '\n';
    char rc;

    while (Serial.available() > 0 && newData == false) {
        rc = Serial.read();
        if (rc != endMarker) {
            if (ndx < numChars - 1) {
                receivedChars[ndx] = rc;
                ndx++;
            }
        } else {
            receivedChars[ndx] = '\0';
            ndx = 0;
            newData = true;
        }
    }
}

// ---------------------------------------------------------------------
// Calibration mode: an open-loop step for identify.cpp.
//
// Fan held OFF throughout, because PLANT_R and PLANT_C describe the
// unforced plant; identifying with the fan running would fold the
// actuator into the model it is supposed to act on.
// ---------------------------------------------------------------------
boolean calMode = false;
unsigned long calStart = 0;
const unsigned long CAL_MS = 700000UL;   // ~700 s, past 5 slow tau

void startCal() {
    calMode = true;
    calStart = millis();
    setFan(0.0);
    setHeater(1.0);
    Serial.println("# thermloop step response");
    Serial.print("# P_watts=");
    Serial.println(P_MAX, 3);
    Serial.println("t,T");
}

// Accepts either the bare [0,1] flag or a "focused"/"relaxed" word, since
// GUI versions have differed on which they emit.
void handleMessage() {
    if (!newData) return;

    if (receivedChars[0] == 'C' && receivedChars[1] == 'A') {
        startCal();
        newData = false;
        return;
    }

    float flag = 0.0;
    if (receivedChars[0] == '1' || receivedChars[0] == 'f'
                                || receivedChars[0] == 'F') {
        flag = 1.0;
    }

    // EMA smoothing. Focus flips fast; the plant does not.
    focusLevel = 0.8 * focusLevel + 0.2 * flag;
    plan.burstNow = 0.30 + 1.30 * focusLevel;

    newData = false;
}

// ---------------------------------------------------------------------
unsigned long lastTick = 0;

void setup() {
    Serial.begin(BAUD);
    pinMode(PIN_HEAT, OUTPUT);
    setHeater(0.0);
    setupFanPWM();

    throttle.begin(PLANT_R, PLANT_C, 2, T_AMB, DT);
    throttle.setPlanStride(STRIDE);
    throttle.setAlphaMin(0.1);
    throttle.setSyncTau(SYNC_TAU);

    Serial.println("<thermloop ready>");
    lastTick = millis();
}

void loop() {
    recvWithEndMarker();
    handleMessage();

    unsigned long now = millis();
    if (now - lastTick < (unsigned long)(DT * 1000.0)) return;
    lastTick = now;

    float T_meas = readTemperature();

    // --- calibration: log the raw step and do nothing else ------------
    if (calMode) {
        unsigned long elapsed = now - calStart;
        if (elapsed > CAL_MS) {
            calMode = false;
            setHeater(0.0);
            Serial.println("# done");
            throttle.reset();
            fan.reset();
            return;
        }
        Serial.print((float)elapsed / 1000.0, 2);
        Serial.print(',');
        Serial.println(T_meas, 3);
        return;
    }

    // --- 1. predict, then throttle the requested workload --------------
    float alpha = throttle.solve(&planAt, &plan, HORIZON, T_LIMIT);
    float P_exec = alpha * planAt(0, &plan);

    // --- 2. drive the heater ------------------------------------------
    setHeater(P_exec / P_MAX);

    // --- 3. advance the model, then re-anchor it on the sensor ---------
    throttle.push(P_exec);
    throttle.syncTemperature(T_meas);

    // --- 4. fan loop closes on the MEASURED temperature ----------------
    //     Error is positive when too hot, so the fan spins up. The
    //     controller clamps to [0,1] with conditional integration, so
    //     sitting at full fan does not wind the integrator up.
    float u = fan.update(T_meas - T_TARGET);
    setFan(u);

    // --- 5. telemetry back to the host --------------------------------
    Serial.print(focusLevel, 2);   Serial.print(',');
    Serial.print(alpha, 3);        Serial.print(',');
    Serial.print(P_exec, 3);       Serial.print(',');
    Serial.print(T_meas, 2);       Serial.print(',');
    Serial.print(throttle.temperature(), 2); Serial.print(',');
    Serial.print(throttle.bias(), 3); Serial.print(',');
    Serial.println(u, 3);
}
