# Feedback Value Certification

Core algorithms for finite-horizon value evaluation of contractive learners whose decisions affect subsequent feedback.

The project provides:

- **Boundary certificates:** group equal projected increments, verify rational word separation, and bound finite-state block survival.
- **Exact scalar values:** propagate reachable envelopes forward and rational value partitions backward, retaining the actual tag at each threshold.
- **RPS value enclosures:** combine an exact prefix with causal grid propagation and directed integer rounding in a two-dimensional rational Rock–paper–scissors/Win-stay–lose-shift model.

[中文说明](README.zh-CN.md)

## Requirements

Python 3.10 or later. The Python implementation uses the standard library and has no runtime package dependencies. Building the RPS executable requires CMake 3.16 or later and a C++17 compiler with `__int128` support, such as GCC or Clang.

## Install

```bash
git clone https://github.com/wyyy5111/feedback-value-certification.git
cd feedback-value-certification
python -m pip install .
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
```

## Run

Check the grouped boundary certificate of the supplied rational RPS kernel:

```bash
feedback-value certificate --input examples/rps_certificate.json
```

The result includes block length `1` and survival upper bound `179/180`, with the projected collision groups and separation data.

Compute an exact scalar expected reward:

```bash
feedback-value scalar --input examples/reactive_binary.json --policy tracking --horizon 6
```

The exact value is `7467729/9765625`. Select the `contrarian` policy to evaluate its own feedback law. `python -m feedback_value_certification` provides the same command-line interface.

Compute deterministic RPS enclosures:

```bash
./build/rps-value --T 6 --G 8
```

The executable writes JSON to standard output and progress messages to standard error. Each endpoint is represented by an integer numerator and a positive denominator. `BR_minus_mode` encloses the difference of the two expectations. The exact prefix has `min(T, 8)` rounds; continuation uses one shared geometric child before the fresh finite-state draw. The fixed rounding denominator is `2^40`.

## Python API

```python
import json
from feedback_value_certification import Model, evaluate_exact

with open("examples/reactive_binary.json", encoding="utf-8") as stream:
    model = Model(json.load(stream), "tracking")
answer = evaluate_exact(model, 6)
print(answer["value"])
```

Model inputs use integers or rational strings such as `"2/3"`. Threshold point tags and open-region tags are specified separately. The response and reward tables define the joint dependence of current reward and the next physical state. Scalar rewards have absolute value at most one.

`evaluate_exact` returns a rational value and construction counts. Its `max_cuts` and `max_atoms` limits default to 12,000 cuts per state/layer and 200,000 constructed atoms. Reaching a limit raises `EvaluationLimitError`; the CLI returns status `RESOURCE_LIMIT` and exit code `2`. Invalid inputs return exit code `1`. A completed certificate query can report an inconclusive sufficient check without an execution error.

The C++ evaluator specializes the rational RPS model with contraction `49/50`, preferred-action probability `58/60`, first-index tie resolution, initial learner point `(1/3, 1/3, 1/3)`, and a uniform initial opponent action. Its memory preflight checks a conservative one-GiB allocation bound. Expected values and intervals refer to this declared rational model.

## Check the implementation

```bash
RPS_VALUE_EXECUTABLE="$PWD/build/rps-value" python -m unittest discover -s tests -v
```

The checks cover exact short-horizon values, threshold tags, rational probability validation, construction limits, causal continuation, and RPS endpoint regressions. GitHub Actions runs the checks on Python 3.10 and 3.12.

## Layout

```text
src/feedback_value_certification/  Python certificate and scalar algorithms
cpp/rps_value.cpp                 C++ RPS enclosure algorithm
examples/                        Two minimal input examples
tests/                           Unit and command-line checks
```

## Source and license

The implementation is available under the [MIT license](LICENSE). [PROVENANCE.md](PROVENANCE.md) records the core source lineage and the upstream mathematical model references. Stable source snapshots are available through the repository's version tags.
