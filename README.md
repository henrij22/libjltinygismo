# libjltinygismo

The C++ component of [TinyGismo.jl](https://github.com/henrij22/TinyGismo.jl) — a
[CxxWrap](https://github.com/JuliaInterop/CxxWrap.jl)-based binding layer over
[G+Smo](https://github.com/gismo/gismo) (Geometry + Simulation Modules), the isogeometric
analysis library.

Building this repository produces a single shared library, `libjltinygismo`, which Julia
loads with `@wrapmodule`. Prebuilt binaries for all major platforms are shipped through
[TinyGismo_jll.jl](https://github.com/henrij22/TinyGismo_jll.jl); you only need to build here
when working on the bindings themselves.

## What it does

The library registers a curated subset of G+Smo with the Julia runtime: knot vectors,
B-spline and NURBS bases (univariate and tensor product), the matching geometries, the
geometry factories, Paraview export and XML file reading.

It is a translation layer, not a thin passthrough. The conventions it establishes are what
make the Julia side feel native:

- **1-based indexing** everywhere it is observable — basis function indices, control point
  indices, knot span indices and parametric directions are shifted on the way in and out.
- **Parametric directions** are 1-based, with `0` meaning "all directions". G+Smo uses
  0-based directions with `-1` for "all", and indexes without bounds checks, so the
  translation happens in one place (`toGismoDir` in `helper.hh`) which rejects out-of-range
  values rather than letting them corrupt memory.
- **Bang methods** (`eval!`, `deriv!`, `active!`, …) replace G+Smo's `_into()` naming and
  take a preallocated `gsMatrix`/`gsVector` as the last argument.
- **Returned arrays are copies.** `toMatrix`, `toVector` and `knotContainer` allocate a
  Julia-owned array rather than viewing the C++ buffer, because the source is almost always
  a temporary that the collector is free to reclaim.
- **Parametric C++ templates map to parametric Julia types**, so `gsTensorBSplineBasis<2>`
  is spelled `TensorBSplineBasis{2}`.

### Layout

| Path | Contents |
| :-- | :-- |
| `src/jltinygismo/gismo.cpp` | module entry point, calls the register functions |
| `src/jltinygismo/helper.hh` | shared conversion helpers between Julia and G+Smo/Eigen types |
| `src/jltinygismo/basis/` | `gsBasis`, B-spline and NURBS bases, knot vectors |
| `src/jltinygismo/geometry/` | `gsGeometry`, B-spline/NURBS curves, tensor surfaces and volumes |
| `src/jltinygismo/utility/` | matrix and vector types, geometry factories, Paraview, file reading |
| `test/` | C++ unit tests (doctest + CTest) |
| `julia/` | Julia integration suite run against the freshly built library |

## Building

### Prerequisites

- A C++20 compiler and CMake ≥ 3.18
- **G+Smo**, built and installed (`v25.07.0` is what CI uses)
- **libcxxwrap-julia**, which ships with the `CxxWrap.jl` package. Ask Julia for its path,
  from an environment that actually has CxxWrap installed — `julia/` in this repository is
  one:

  ```bash
  julia --project=julia -e 'using Pkg; Pkg.instantiate()'
  julia --project=julia -e 'using CxxWrap; println(CxxWrap.prefix_path())'
  ```

Building G+Smo from source, if you do not have it:

```bash
git clone --depth 1 --branch v25.07.0 https://github.com/gismo/gismo.git
cmake -S gismo -B gismo/build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=$HOME/dev/install \
  -DGISMO_BUILD_LIB=ON
cmake --build gismo/build --parallel
cmake --install gismo/build
```

### Configure and build

Point `CMAKE_PREFIX_PATH` at both dependencies:

```bash
CXXWRAP_PREFIX=$(julia --project=julia -e 'using CxxWrap; println(CxxWrap.prefix_path())')

cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$HOME/dev/install;$CXXWRAP_PREFIX"
cmake --build build --parallel
```

The result is `build/lib/libjltinygismo.so` (`.dylib` on macOS, `.dll` on Windows).

### CMake presets

`CMakePresets.json` carries the `debug`, `release` and `relwithdebinfo` presets. They hard-code
local paths, so adjust `CMAKE_PREFIX_PATH` and `CMAKE_INSTALL_PREFIX` in that file for your
machine, then:

```bash
cmake --preset debug
cmake --build --preset debug
```

### Options

| Option | Default | Meaning |
| :-- | :-- | :-- |
| `JLTINYGISMO_BUILD_TESTS` | `ON` | Build the C++ unit tests. Fetches doctest if it is not installed, so turn this off for a network-free build. |

## Testing

There are two suites, covering different things.

### C++ unit tests

```bash
ctest --test-dir build --output-on-failure
```

These cover the parts of `helper.hh` that do not touch the Julia runtime — the direction
translation and the index shifting. Everything else in this repository is CxxWrap
registration code, which cannot be meaningfully exercised without a live Julia runtime, so it
is covered by the suite below instead. Add new cases to `test/test_helper.cpp`, or a new file
listed in `test/CMakeLists.txt`.

### Julia integration tests

This is the real coverage. It loads the library you just built and drives the bindings
exactly as TinyGismo.jl does:

```bash
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
julia --project=julia julia/runtests.jl build/lib
```

The path argument is the directory containing the shared library; alternatively set
`JLTINYGISMO_LIB`. If G+Smo is installed somewhere the loader does not search, prefix the
command with `LD_LIBRARY_PATH=$HOME/dev/install/lib` (`DYLD_LIBRARY_PATH` on macOS).

The suite checks the conventions above (1-based indices, direction handling, array
ownership), verifies mathematical properties that must hold regardless of implementation
(partition of unity, refinement preserving the geometry, NURBS circles having radius exactly
one, derivatives matching finite differences), and carries a regression test for every fixed
binding bug.

## Pointing TinyGismo.jl at a local build

The suite above loads the library by path, but `TinyGismo.jl` normally resolves it through
the released `TinyGismo_jll`. To run TinyGismo.jl — its test suite, its documentation, or
just a REPL session — against a library you built here, override the JLL with a
`LocalPreferences.toml` next to the `Project.toml` of whichever environment you are using:

```toml
[TinyGismo_jll]
libjltinygismo_path = "/absolute/path/to/libjltinygismo/build/lib/libjltinygismo.so"
```

`JLLWrappers` reads that preference instead of the artifact, so nothing else changes:

```bash
cd /path/to/TinyGismo.jl
cat > LocalPreferences.toml <<'EOF'
[TinyGismo_jll]
libjltinygismo_path = "/absolute/path/to/libjltinygismo/build/lib/libjltinygismo.so"
EOF

julia --project=. -e 'using Pkg; Pkg.test()'
```

Notes:

- The path must be the **shared library file itself**, absolute, including the extension
  (`.so`, `.dylib`, `.dll`).
- The preference is read at precompilation time. Julia notices the change and recompiles, but
  if a stale binary seems to persist, delete the package's compile cache or touch a source
  file to force it.
- Each environment needs its own file. Building the docs also means dropping one in `docs/`.
- Add `LocalPreferences.toml` to `.gitignore` in the consuming repository — it holds an
  absolute path from your machine and must never be committed.
- Once the JLL is released with your changes, delete the file so the environment goes back to
  the published binary. A forgotten override is a confusing way to lose an afternoon.

This is also the way to reproduce a bug against an unreleased fix: build here, point
TinyGismo.jl at the result, and run its tests.

## Working on the bindings

A new binding is a `method` call on the relevant type wrapper. Points worth remembering:

- Shift every index that a user can observe. `i - 1` going in, `incrementByOne` coming out.
- A lambda's parameter list and its `arg(...)` specifiers must agree in arity. A mismatch
  compiles silently and registers a Julia method that cannot be called correctly.
- Return `gsMatrix`/`gsVector` by value and let `toMatrix`/`toVector` copy it out. Do not
  hand back `jlcxx::make_julia_array` over memory that C++ owns.
- Types exposed to Julia that inherit from `gsBasis` or `gsGeometry` need a
  `jlcxx::SuperType` specialization, placed **before** the `add_type` call that registers
  them, or the inherited methods fail at runtime with an upcast error.
- Generic lambdas (`[](auto& x, ...)`) cannot be wrapped; spell out the concrete type.

After changing a binding, run both suites, and add a case to `julia/runtests.jl` covering the
new behavior.

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).
