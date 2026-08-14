# Hill four-body Poincaré sections

This repository computes Poincaré sections for the McGehee-regularized,
planar Hill four-body problem with oblateness. The checked-in parameters are a
Sun–Jupiter–Hektor–Skamandrios study case. In practice, the code samples initial
conditions on the section `y = 0`, follows the regularized flow, and records
positive crossings (`ydot > 0`).

There are three ways to use it:

- `hill4bp_poincare` is a small C++17 program with its own floating-point ODE
  integrator. It is the quickest way to build the model, run a smoke test, and
  make an ordinary Poincaré plot.
- `hill4bp_poincare_capd` uses CAPD's Taylor solver for larger floating-point
  surveys and records a detailed outcome for every initial condition.
- `hill4bp_validate_capd` uses CAPD interval arithmetic and C0 sets to validate
  discrete points or covers near the first neck. This is the proof-oriented
  program; the output of the other two programs is numerical evidence, not a
  computer-assisted proof.

The source, build files, and working directory for all three programs are in
[`Hill_4BP_cpp_1`](Hill_4BP_cpp_1/). The sibling [`CAPD`](CAPD/) directory is an
upstream checkout of CAPD::DynSys, not part of the Hill 4BP implementation.

![Normalized first-return sections near the E1 neck](Hill_4BP_cpp_1/professional_neck_visualizations/figures/03_normalized_first_returns.png)

## Quick start

For a first run, only a C++17 compiler and either Make or CMake are needed.
Run these commands from the top of this repository:

```sh
cd Hill_4BP_cpp_1
make -B
./build/hill4bp_poincare --quick --quiet \
  --output quick_section.csv \
  --svg quick_section.svg \
  --metadata quick_section.txt
```

The command writes the full-precision section points to `quick_section.csv`, a
browser-ready preview to `quick_section.svg`, and the resolved parameters and
run statistics to `quick_section.txt`. The `--quick` preset is deliberately
small. Omit it to use the grid in `Hill_4BP_parameters.cfg`; the default grid
is much slower.

To see every command-line option:

```sh
./build/hill4bp_poincare --help
```

## Dependencies

| Task | Required | Optional |
| --- | --- | --- |
| Build and run `hill4bp_poincare` | C++17 compiler; Make or CMake 3.16+ | Bash for the shell tests |
| Build the CAPD programs | The above, plus a built CAPD::DynSys tree and its `capd-config` script | Boost and the other CAPD build prerequisites, depending on the CAPD configuration |
| View CSV results | Python 3, NumPy, Matplotlib | A GUI backend for an interactive window |
| Run the original/reference workflow | MATLAB | None of the C++ programs require MATLAB |
| Regenerate the paper figures | MATLAB | LaTeX for the accompanying manuscript |

The included CAPD tree can be built with native interval arithmetic and
multiprecision disabled. That is the least complicated configuration for this
project. On Ubuntu or WSL, its usual prerequisites can be installed with:

```sh
sudo apt update
sudo apt install -y build-essential cmake git pkg-config autoconf libtool \
  libboost-all-dev
```

On macOS, install the Xcode command-line tools and CMake. For an Apple Silicon
build, explicitly request `arm64`; otherwise a CMake process running under
Rosetta can produce an x86_64 CAPD library. CAPD's bundled FILIB backend is not
used here.

Python plotting support can be installed into an environment of your choice:

```sh
python3 -m pip install numpy matplotlib
```

## Building

### Standalone program

The Makefile is the shortest route:

```sh
cd Hill_4BP_cpp_1
make -B
make test
```

`-B` rebuilds the executable for the current machine. Files already present in
`build/` are saved artifacts and should not be assumed to match the host
architecture or local libraries.

The equivalent CMake build is:

```sh
cd Hill_4BP_cpp_1
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

Both produce `build/hill4bp_poincare`. The tests cover parameter parsing,
coordinate conversions, the regularized vector field, conserved quantities,
equilibria, and the shell-study contracts.

### CAPD survey and validator

First build the included CAPD checkout. From the top of this repository:

```sh
cmake -S CAPD -B CAPD/build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCAPD_INTERVAL_TYPE=NATIVE \
  -DCAPD_ENABLE_MULTIPRECISION=OFF
cmake --build CAPD/build -j
```

On an Apple Silicon Mac, add
`-DCMAKE_OSX_ARCHITECTURES=arm64` to the configure command. You can check the
result with `lipo -archs CAPD/build/libcapd.a`.

Then build both Hill 4BP CAPD programs:

```sh
cd Hill_4BP_cpp_1
export CAPD_CONFIG="$PWD/../CAPD/build/bin/capd-config"
bash apps/build_capd_apps.sh
```

The build script uses the same CAPD installation for both programs, supplies
rounding-safe compiler flags, and performs a dry-run preflight. If CAPD is
installed elsewhere, set `CAPD_CONFIG` to that installation's executable
`capd-config`. `make capd` is a convenience wrapper around the same script.

## Running the programs

All examples below assume the current directory is `Hill_4BP_cpp_1`. Running
from there also ensures that the programs find `Hill_4BP_parameters.cfg`.

### Standalone Poincaré section

```sh
./build/hill4bp_poincare \
  --x-count 120 \
  --xdot-count 21 \
  --iterates 100 \
  --output poincare_grid.csv \
  --svg poincare_section.svg \
  --metadata poincare_metadata.txt
```

This program uses conventional `--option value` arguments. Command-line
values override the shared configuration for that run only. Use `--no-svg`
when the point cloud is large, or `--svg-stride N` to put only every Nth point
in the preview while retaining every point in the CSV.

### CAPD floating-point survey

CAPD programs use case-sensitive `key=value` arguments, with no leading
dashes. A dry run is a useful first step because it resolves the equilibrium
energy, neck window, and grid without integrating or creating output files:

```sh
./build/hill4bp_poincare_capd \
  neck=E1 side=open delta=1e-4 \
  neck_window_sigma=4 \
  x_count=60 xdot_count=11 iterates=5 \
  dry_run=1
```

Remove `dry_run=1` and name the outputs to run the survey:

```sh
./build/hill4bp_poincare_capd \
  neck=E1 side=open delta=1e-4 \
  neck_window_sigma=4 \
  velocity_mode=fraction \
  x_count=60 xdot_lower=-0.9 xdot_upper=0.9 xdot_count=11 \
  iterates=5 \
  output=e1_open_points.csv \
  returns=e1_open_returns.csv \
  outcomes=e1_open_outcomes.csv \
  zvb=e1_open_zvb.csv \
  svg=e1_open.svg \
  metadata=e1_open.txt
```

Here `delta` is the positive Hamiltonian distance from the selected critical
neck. `side=closed` uses `h = H(E1) - delta`; `side=open` uses
`h = H(E1) + delta`. The Jacobi constant is `C = -2h`. Do not run exactly at
`delta=0`: the equilibrium has zero speed and is non-transverse to this
section.

The two-sided reconnaissance script runs both sides for several values of
`delta` and keeps separate CSV, SVG, metadata, outcome, and log files:

```sh
DELTA_LEVELS="1e-2 1e-3 1e-4" \
X_COUNT=60 XDOT_COUNT=11 ITERATES=5 \
bash apps/run_neck_sweep_capd.sh
```

### Rigorous CAPD validation

The interval validator currently supports the E1/E2 neck. Start with its
preflight:

```sh
./build/hill4bp_validate_capd \
  neck=E1 side=open delta=1e-4 \
  neck_window_sigma=4 \
  grid_mode=points x_count=21 fraction_count=11 \
  iterates=1 dry_run=1
```

For actual proof studies, use `apps/run_professional_neck_study.sh`. It splits
the grid into deterministic shards, records configuration and binary hashes,
keeps the per-shard logs, and merges results only after all shards finish. A
small pilot looks like this:

```sh
OUTPUT_DIR="$PWD/professional_neck_studies/pilot" \
MODE=both \
DELTA_LEVELS="1e-3 1e-4" \
SURVEY_X_COUNT=40 SURVEY_FRACTION_COUNT=15 SURVEY_ITERATES=5 \
PROOF_GRID_MODE=points PROOF_X_COUNT=11 PROOF_FRACTION_COUNT=7 \
PROOF_ITERATES=1 PROOF_REQUIRE_RETURNS=0 REQUIRE_COMPLETE=1 \
SHARD_COUNT=4 JOBS=4 \
bash apps/run_professional_neck_study.sh
```

Interval validation is substantially more expensive than a floating-point
survey. Read the [detailed project notes](Hill_4BP_cpp_1/README.md) before
turning a pilot into a theorem claim; in particular, an unresolved leaf is a
remaining proof obligation. A point grid certifies the sampled points, not the
rectangles between them. Use `grid_mode=cover` only for a carefully chosen,
strictly admissible parameter rectangle.

### Plotting a CSV

```sh
python3 view_poincare.py e1_open_points.csv
python3 view_poincare.py e1_open_points.csv --png e1_open.png
python3 view_poincare.py e1_open_points.csv --pdf e1_open.pdf
```

These examples use the headerless, two-column point cloud written by the CAPD
survey. The viewer keeps the source coordinates in memory and redraws them
while zooming. The generated SVG is a convenient preview; the CSV is the
numerical result to keep.

## Changing constants and defaults

Edit [`Hill_4BP_parameters.cfg`](Hill_4BP_cpp_1/Hill_4BP_parameters.cfg). The
C++ and CAPD programs read it at startup, so changing the file does not require
a rebuild.

The file has two groups of values:

- `G`, `mu`, `u1`, `u2`, `c1`, `c2`, `c3`, `lambda1`, `lambda2`, `nu`, and
  `alpha` define the model. If `lambda1` or `lambda2` is omitted, the standalone
  library and floating-point CAPD survey derive it from `mu`, `u1`, and `u2`;
  the rigorous validator requires both decimal values to be present.
- Names beginning with `poincare_` set the default energy, grid bounds and
  counts, return count, section offset, solver tolerances, step sizes, and SVG
  stride.

Most run settings can also be overridden on the command line. The standalone
program uses names such as `--x-count 120`; the CAPD survey uses names such as
`x_count=120`. Run either program with `--help` for the full list.

There is one important model boundary in the present CAPD implementation. Its
regularized field is specialized to:

```text
c1 = 0
c2 = 0
alpha = 3
nu = 1
```

The checked-in `c3`, `lambda1`, and `lambda2` are intentional continuation
inputs. Do not set `c1` or `c2` to a nonzero value and assume that all three
bodies have thereby been modeled as independently oblate. The Hill-limit
coefficients and CAPD field must first be extended; the CAPD programs reject
unsupported values rather than silently run a different model.

For a one-off energy or grid change, prefer a command-line override. For a
reproducible change to the physical system or shared defaults, edit the config
and keep the metadata written by the run.

## Repository map

```text
.
├── README.md                         this file
├── CAPD/                            upstream CAPD::DynSys checkout
└── Hill_4BP_cpp_1/
    ├── apps/                        executables and batch scripts
    ├── include/crtbp_cpp/           public C++ headers
    ├── src/                         standalone model and ODE code
    ├── tests/                       C++ regressions and shell tests
    ├── professional_neck_studies/   saved survey/proof products
    ├── professional_neck_visualizations/
    ├── *.m                          MATLAB reference implementation
    ├── Hill_4BP_parameters.cfg      model and run defaults
    ├── CMakeLists.txt / Makefile
    └── view_poincare.py
```

The [visualization notes](Hill_4BP_cpp_1/professional_neck_visualizations/README.md)
explain the saved E1 study, normalization, figure meanings, and the limits of
what those results establish.

## References and citation

The regularized vector field follows equation (4.9) of Belbruno, Gidea, and
Lam, [*Regularization of the Hill four-body problem with oblate
bodies*](https://doi.org/10.1007/s10569-023-10122-x), *Celestial Mechanics and
Dynamical Astronomy* **135**, 6 (2023). The physical study case is described by
Burgos-García et al., [*Hill Four-Body Problem with Oblate Bodies: An
Application to the Sun–Jupiter–Hektor–Skamandrios
System*](https://doi.org/10.1007/s00332-020-09640-x), *Journal of Nonlinear
Science* **30**, 2925–2970 (2020).

The CAPD-backed programs depend on [CAPD::DynSys](https://github.com/CAPDGroup/CAPD).
If results from either CAPD program are used in a paper, cite:

> Tomasz Kapela, Marian Mrozek, Daniel Wilczak, and Piotr Zgliczyński,
> “CAPD::DynSys: A flexible C++ toolbox for rigorous numerical analysis of
> dynamical systems,” *Communications in Nonlinear Science and Numerical
> Simulation* **101** (2021), 105578.
> [doi:10.1016/j.cnsns.2020.105578](https://doi.org/10.1016/j.cnsns.2020.105578)

```bibtex
@article{Kapela2021CAPDDynSys,
  author  = {Kapela, Tomasz and Mrozek, Marian and Wilczak, Daniel and
             Zgliczy{\'n}ski, Piotr},
  title   = {{CAPD::DynSys}: A Flexible {C++} Toolbox for Rigorous Numerical
             Analysis of Dynamical Systems},
  journal = {Communications in Nonlinear Science and Numerical Simulation},
  volume  = {101},
  pages   = {105578},
  year    = {2021},
  doi     = {10.1016/j.cnsns.2020.105578}
}
```

## License

The Hill 4BP code in `Hill_4BP_cpp_1` is released under the
[MIT License](Hill_4BP_cpp_1/LICENSE). CAPD is maintained and licensed
separately by the CAPD Group; see [`CAPD/COPYING`](CAPD/COPYING) and the
copyright notices in that source tree.
