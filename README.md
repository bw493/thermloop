# thermloop

A C++ thermal modeling and control platform built on Foster RC networks and the
electro-thermal analogy. thermloop predicts a device's temperature trajectory
from its planned workload, throttles that workload *before* a thermal limit is
crossed, and closes a PI fan loop on the measured temperature — all from a
plant model fitted to your actual hardware with a single step response.

The intended application is thermal management of wearable and embedded BCI
acquisition hardware, where sustained DSP bursts (filtering, classification)
heat small enclosures against skin-contact temperature limits. The onboard
core is freestanding — no heap, no exceptions, no STL — and runs on anything
from an ATmega328P upward.

## Supported targets

| Target | Precision | Max RC stages | Status |
|---|---|---|---|
| Arduino UNO Rev3 (ATmega328P) | `float` | 4 | Complete reference sketch: `focus_throttle.ino` |
| ESP32-S3 (ADS1299-based EEG boards) | `float` (LX7 FPU is single-precision) | 8 | Portable core; adapt the reference sketch (see below) |
| Raspberry Pi / Linux ARM (PiEEG) | `double` | 16 | Portable core; build as ordinary Linux C++ |

The platform profile is selected automatically in `throttle.h` from the
compiler's architecture defines. To override, define `TL_REAL` and
`TL_MAX_STAGES` before including the header.

## Repository map

**Onboard core (portable, freestanding)**
- `throttle.h` — `ThermalThrottle`: exact recursive state-space realization of
  the Foster plant (one scalar per stage instead of a multi-kilobyte
  convolution kernel), predictive `solve()` returning a throttle factor
  `alpha`, a planning-stride knob that cuts compute cost without degrading the
  physics, an integral disturbance estimator (`syncTemperature()`) that
  re-anchors the model on a real sensor, and `solveBlocks()` for classifier
  pipelines that drop whole DSP blocks rather than scaling a clock.
- `controller.h` — Tustin-discretized `PIController` with output clamping and
  conditional (anti-windup) integration, plus `ContinuousPIReference`, a
  fine-grained continuous integrator used only as validation ground truth.

**Host-side modeling and design tools (desktop C++)**
- `plant.h` — `FosterNetwork`: the shared plant model (impulse, step, time
  constants). Every other module composes against this header.
- `identify.cpp` — fits `PLANT_R` / `PLANT_C` from a measured step response by
  variable-projection least squares: linear in R, so only the time constants
  are grid-searched. Includes `--selftest` and an NTC beta two-point solver.
- `fanplant.h` — fan actuator model: airflow scales the convection (last
  Foster) stage, `R_conv(u) = R_conv0 / (1 + kFan * u)`.
- `loopshape.h` — `LoopShaper`: crossover frequency and phase margin of
  `L(s) = C(s) * H_fan(s)`, and `suggestPI()`, which inverts the problem —
  choose a crossover and margin, receive the `(Kp, Ki)` that deliver them.
- `impedance.h` — `ThermalImpedance`: the driving-point impedance
  `Z(s) = sum_i R_i / (1 + s * tau_i)` evaluated at complex `s`.
- `predict.h` / `predict.cpp` — convolution-based temperature predictor
  (the kernel-storage realization; `throttle.h` is the memory-lean equivalent).

**Reference integration**
- `focus_throttle.ino` — Arduino UNO Rev3 sketch closing the platform around
  the OpenBCI GUI Focus Widget over serial.

**Tests and pedagogy**
- `throttle_test.cpp`, `loopshape_test.cpp`, `controller_test.cpp` — desktop
  validation drivers (equivalence of the two plant realizations, stride
  exactness, sync convergence, discrete-vs-continuous PI agreement).
- `plant.cpp`, `plant_single.cpp`, `advisor.cpp` — the teaching track: a
  standalone plant demo, a trapezoidal-integration walkthrough, and the
  coffee-cooling advisor that introduces the plant/controller composition.

## Quick start (any desktop)

```sh
g++ -std=c++17 -O2 throttle_test.cpp   -o throttle_test   && ./throttle_test
g++ -std=c++17 -O2 loopshape_test.cpp  -o loopshape_test  && ./loopshape_test
g++ -std=c++17 -O2 controller_test.cpp -o controller_test && ./controller_test
g++ -std=c++17 -O2 identify.cpp        -o identify        && ./identify --selftest
```

All four should pass before you trust anything on hardware.

## Workflow: from your rig to a running loop

The constants shipped in `focus_throttle.ino` are self-test defaults and
**will be wrong for your hardware**. The platform is designed so you never
have to guess them:

1. **Calibrate.** Flash the sketch, send `CAL` over serial, and capture the
   CSV it streams (an open-loop heater step, fan held off, roughly 700 s).
2. **Identify.** `./identify cal.csv <P_watts>` fits the two-stage Foster
   model and prints the `PLANT_R` / `PLANT_C` lines to paste back into the
   sketch (or into your own firmware).
3. **Design the fan loop.** Point `LoopShaper::suggestPI()` at the fitted
   plant via `fanplant.h`, choose a crossover and phase margin, and hand the
   resulting `(Kp, Ki)` to `PIController`. `controller_test.cpp` shows the
   Tustin discretization agreeing with the continuous reference.
4. **Deploy.** `ThermalThrottle::begin()` with the fitted R, C at your control
   interval; call `solve()` each tick against your workload plan; `push()` the
   executed power; `syncTemperature()` with the measured value.

## Target-specific notes

### Arduino UNO Rev3 + OpenBCI GUI Focus Widget

`focus_throttle.ino` is complete and wired as follows:

- **A0** — NTC divider (10 k NTC to GND, 10 k fixed to 5 V)
- **D6** — MOSFET gate driving the power resistor (the "workload")
- **D9** — fan PWM, reconfigured to 25 kHz via Timer1 (`ICR1` top; note this
  removes D9/D10 from `analogWrite()` — `millis()` on Timer0 is unaffected)

Serial protocol follows the OpenBCI reference sketch: **57600 baud,
newline-terminated, 32-byte buffer**. The sketch accepts either the Focus
Widget's bare `0`/`1` flag or the `focused`/`relaxed` words, EMA-smooths the
flag into a continuous demand level, and re-solves the predictive throttle
every 0.5 s over a 60 s lookahead. The host performs all EEG DSP; the UNO
never sees a sample of EEG — focus arrives purely as a workload demand.

### ESP32-S3 (ADS1299-based boards)

`throttle.h` and `controller.h` compile unmodified under the Arduino-ESP32
core; the platform profile selects single-precision `float` (the Xtensa LX7
FPU emulates `double` in software) and raises the stage ceiling to 8. Two
adaptations are required when porting the reference sketch:

- The 25 kHz fan PWM in `setupFanPWM()` writes AVR Timer1 registers directly
  and will not compile on Xtensa. Replace it with the ESP32 LEDC peripheral
  at 25 kHz, and replace `analogRead()` scaling with the S3's 12-bit ADC.
- If your acquisition pipeline degrades by dropping whole DSP blocks (e.g.
  disabling classifier stages) rather than scaling a clock, call
  `solveBlocks(plan, ctx, K, T_limit, nBlocks)` instead of `solve()`; it
  quantizes the throttle onto a `k/nBlocks` grid, rounding down so the
  temperature limit remains a ceiling.

### Raspberry Pi (PiEEG)

On Linux/ARM the core builds as ordinary C++ at full `double` precision with
up to 16 stages — no Arduino toolchain involved:

```cpp
#include "throttle.h"

ThermalThrottle throttle;

int main() {
    const double R[] = {4.15, 8.88};   // from identify.cpp, not from a guess
    const double C[] = {1.59, 14.61};
    throttle.begin(R, C, 2, 24.0, 0.5);
    throttle.setSyncTau(600.0);

    // each control tick:
    //   double alpha = throttle.solve(plan_fn, &ctx, horizon, T_limit);
    //   ...execute alpha-scaled workload...
    //   throttle.push(P_executed);
    //   throttle.syncTemperature(read_sensor());
}
```

Read the temperature from whatever sensor the rig exposes (an I2C sensor or
the SoC thermal zone under `/sys/class/thermal/`), and run the calibration and
identification workflow above exactly as on the microcontrollers — the model
is only as good as the step response behind it.

## Design principles

- **One system, two realizations.** The convolution predictor (`predict.h`)
  and the recursive state-space throttle (`throttle.h`) are the same LTI
  system exactly; only the memory bill differs. `throttle_test.cpp` verifies
  their equivalence term by term.
- **Measure, then model.** Every plant constant traces to a fitted step
  response; `identify.cpp` exists so no parameter is a plausible-looking guess.
- **Ceilings, not targets.** The throttle rounds conservatively and plans
  against a rising worst-case demand envelope, because over-throttling costs
  a little responsiveness while under-throttling costs an overshoot.

## License

MIT — see `LICENSE`.
