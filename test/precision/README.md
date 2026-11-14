# GLM precision regression

This suite pins numerical behavior for view/projection matrices, matrix
decomposition, Euler-angle round trips, quaternion matrix casts and `slerp`.
It runs the same checks for `float` and `double`.
The implementation is `test/gtx/gtx_precision_regression.cpp`.

Run it from any directory with:

```bash
test/precision/run.sh
```

The script configures `GLM_BUILD_TESTS=ON`, builds only
`test-gtx_precision_regression`, and runs the matching CTest case. Set
`GLM_PRECISION_BUILD_DIR` to use a different build directory.

Equivalent manual commands:

```bash
cmake -S . -B build-precision -DGLM_BUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-precision --target test-gtx_precision_regression --parallel
ctest --test-dir build-precision -R test-gtx_precision_regression --output-on-failure
```

## Fixed inputs and golden values

Projection inputs are `fovy=1.2`, `aspect=1.5`, `zNear=0.25`, `zFar=100`.
The shared elements are:

| Element | Value |
| --- | --- |
| `m[0][0]` | `0.9744639647187348` |
| `m[1][1]` | `1.4616959470781021` |
| NO `m[3][2]` | `-0.5012531328320802` |
| ZO `m[3][2]` | `-0.2506265664160401` |

Depth elements are `-1.005012531328321` for RH_NO, `-1.0025062656641603`
for RH_ZO, `1.005012531328321` for LH_NO and `1.0025062656641603` for
LH_ZO. The `m[2][3]` sign is RH `-1`, LH `+1`.

The extreme-depth cases use `zNear=1e-6` and `zFar=numeric_limits<T>::max()`.
Their depth translation must remain finite at `-2e-6` for NO and `-1e-6`
for ZO. This specifically catches multiplication orders such as
`2 * zFar * zNear` overflowing before the division.

View inputs are `eye=(1.3,-2.2,4.7)`, `center=(0.4,1.1,0.2)` and
`up=(0,1,0)`. RH must map `center` to `(0,0,-5.6524331044250316)`, LH to
`(0,0,+5.6524331044250316)`; both 3x3 rotation blocks must be orthonormal
with determinant `+1`.

The near-half-turn slerp case uses the unit axis
`normalize(vec3(1, 2, 3))` and angle `pi - 1e-4`. At `a=0.5` the golden
quaternion is:

```text
w=0.7071244586351045
x=0.18897751188964454
y=0.3779550237792891
z=0.5669325356689336
```

Matrix decomposition uses a fixed translation, `YXZ` rotation, shear and
negative scale `(-2, 0.5, 3)`. The suite checks recomposition, the scale/determinant
product and reflection sign rather than requiring a unique sign assignment for
negative scale components.

## Precision thresholds

Errors are normalized as
`abs(actual - expected) / (absolute + relative * abs(expected))`; a passing
check has normalized error at most `1`.

| Check | float absolute | float relative | double absolute | double relative |
| --- | ---: | ---: | ---: | ---: |
| Projection elements | `2e-7` | `2e-6` | `2e-14` | `2e-14` |
| Extreme projection depth | `1e-10` | `2e-5` | `1e-18` | `2e-14` |
| View orthogonality | `5e-6` | - | `2e-14` | - |
| View/projection determinant | `5e-6` | `5e-6` | `2e-14` | `2e-14` |
| Decompose/recompose | `3e-6` | `3e-6` | `2e-14` | `2e-14` |
| Small-scale decompose (`1e-8`) | `1e-12` | `3e-6` | `1e-20` | `2e-14` |
| Euler angle/matrix round trip | `2e-6` | `2e-6` | `2e-14` | `2e-14` |
| Quaternion component/angle | `6e-6` | `6e-6` | `3e-14` | `3e-14` |

Not-a-number and infinity are handled as boundary outcomes, not tolerances:

- `decompose` rejects any NaN or infinity in the input and returns `false`.
- Projection and `slerp` inputs propagate NaN/infinity into non-finite output;
  the tests assert non-finiteness rather than accepting an arbitrary finite value.
- `slerp` normalizes finite non-unit inputs before interpolation, so the result
  is unit length and independent of positive input scale.

The executable prints both threshold sets before running. Each failure reports
the test name, metric, measured normalized error and limit.
