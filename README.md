# sl3d

This repository contains codes for forward and adjoint modelling of glacial isostatic adjustment (GIA). Programs contained in gia1d are for simulations in 1D viscosity models, but with the ability to compute sensitivity kernels for 3D viscosity perturbations [[1]](#1).
Programs in sl3d extend the methods to allow for 3D viscosity models using a pseudospectral method described in [[2]](#2). These codes were originally developed by Ophelia Crawford and David Al-Attar, and have been subsequently extended as part of the
grant NE/V010433/1:NSFGEO-NERC: Adjoint tomography of mantle viscosity using deglacial sea level observations.

The two codes are independent of one another: they share much of their source
history but have since diverged, and each builds its own copy of the modules.


## Repository layout

| Path | Contents |
| --- | --- |
| `gia1d/src` | Sources for the 1D viscosity code |
| `sl3d/src`  | Sources for the 3D viscosity code |
| `data` | Input data: the example earth model, S20RTS, coastlines, and the ice histories |
| `scripts` | Python scripts for plotting the results |
| `build/run` | Created by the build: the place to run the codes from |
| `build/lib` | Created by the build: the two core libraries |

Within each `src` directory the physics lives in a set of modules that are
compiled into a static library, and each program is a thin driver on top of it:

| Module | Role |
| --- | --- |
| `module_model` | Earth model: reading, splining, spectral element meshing, 3D viscosity |
| `module_fourier` | Spherical harmonic transforms, grids and rotations |
| `module_mat` | Assembly of the spectral element system matrices |
| `module_gia` | Time and load interpolation, surface fields, map projections, and in sl3d the Runge-Kutta update of the viscoelastic response |
| `module_sl` | Sea level equation, ice loading, and the forward solver `forward_full_sl` |
| `module_ice` | The choice of ice model, and reading the ice histories from netCDF |
| `module_kern` | Adjoint calculations, sensitivity kernels and inversions |
| `module_util`, `module_spline`, `nrtype`, `nrutil` | Numerical and I/O utilities |


## Building

The build needs gfortran, CMake (>= 3.14), LAPACK/BLAS, FFTW3, OpenMP and
netCDF-Fortran. On Debian or Ubuntu:

```
sudo apt install gfortran cmake liblapack-dev libfftw3-dev libnetcdff-dev
```

Both codes are then built together from the top of the repository:

```
cmake -S . -B build
cmake --build build -j
```

This creates `build/run`, which is where the codes are meant to be run from.
It holds the executables, in a directory per code, beside the input data:

```
build/run/
├── gia1d/          sea_level_1d, ice_kern_1d, visc_kern_1d, ...
├── sl3d/           sea_level_3d, ice_kern_3d, visc_kern_3d, ...
├── prem.simp.vm1   the example earth model, copied from data
├── S20RTS.sph      the shear wave model sl3d scales viscosity from
├── coast_lines_110_rob.dat
└── ice5g, ice6g, ice7g -> links to the ice histories in data
```

The two codes share a few program names, which is why each keeps its own
directory. Nothing is written back into the source tree, so deleting `build`
removes every build product and every result. To put the executables somewhere
permanent instead:

```
cmake --install build --prefix ~/.local
```

Useful options:

| Option | Default | Effect |
| --- | --- | --- |
| `-DCMAKE_BUILD_TYPE=` | `Release` | Standard CMake build types |
| `-DSL3D_BOUNDS_CHECK=OFF` | `ON` | Turns off runtime array bounds checking |

Bounds checking is on by default because the original Makefiles always compiled
with `-fbounds-check`. It costs a fair amount of run time, so it is worth
switching off once a calculation is known to be behaving.


## Programs

### Running them

The programs take the earth model as their first argument, and optionally the
ice model as a second. Everything else they read, and everything they write,
goes in the current working directory, with no path arguments anywhere. So the
place you run from matters, and `build/run` is that place:

```
cd build/run
./gia1d/sea_level_1d prem.simp.vm1              # 1D viscosity, ICE5G
./sl3d/sea_level_3d prem.simp.vm1 ice6g         # 3D viscosity, ICE6G
```

Run with no arguments to see what a program expects. Results land beside the
executables in `build/run`, which is inside the build tree, so they are never
mixed into the repository and `rm -rf build` clears them.

`build/run` is only a convenience, not a requirement: the programs will run in
any directory that has the inputs they need. That means the model file you name
on the command line, `S20RTS.sph` for the sl3d programs that scale viscosity
from it, and an `ice5g`, `ice6g` or `ice7g` directory of ice histories, unless
`ICE_DIR` names the directory that holds them. For a run you want to keep separate:

```
mkdir -p ~/runs/vm1 && cd ~/runs/vm1
cp ~/dev/sl3d/data/prem.simp.vm1 ~/dev/sl3d/data/S20RTS.sph .
ICE_DIR=~/dev/sl3d/data ~/dev/sl3d/build/run/sl3d/sea_level_3d prem.simp.vm1 ice6g
```

They are OpenMP parallel, so set `OMP_NUM_THREADS` to suit the machine.
Results do not depend on how many threads are used: the assembly loops
take the elements in two passes, odd numbered and then even, so that no
two elements sharing a node are ever handled at the same time.

The forward programs write one file per output time, named for the age it
holds: `sl.1d.20.50` from gia1d, `sl.3d.20.50` from sl3d, both meaning 20.5 ka.
The age comes from the run's own time base, so it follows the ice model and
however often the run writes output.


### gia1d — 1D viscosity models

| Program | Purpose |
| --- | --- |
| `sea_level_1d` | Forward calculation of the sea level response to a deglaciation history |
| `visc_kern_1d` | Adjoint sensitivity kernel of a point sea level measurement with respect to 3D viscosity |
| `ice_kern_1d` | Adjoint sensitivity kernel of a point sea level measurement with respect to ice thickness |
| `grace_visc_kern` | As `visc_kern_1d`, but for a GRACE-like spherical harmonic coefficient of the gravity field |
| `sl_init_inv` | Inversion for the initial sea level, given final sea level data |
| `merge_files` | Combines two output fields onto a common grid for plotting |

### sl3d — 3D viscosity models

| Program | Purpose |
| --- | --- |
| `sea_level_3d` | Forward calculation in a laterally varying viscosity model |
| `visc_kern_3d` | Viscosity sensitivity kernel for a point sea level measurement |
| `ice_kern_3d` | Ice thickness sensitivity kernel for a point sea level measurement |
| `grace_visc_inv` | Viscosity inversion driven by GRACE gravity data |
| `sl_visc_inv` | Viscosity inversion driven by sea level data |
| `rsl_inversion` | Viscosity inversion driven by a relative sea level data set |
| `visco_test` | Time steps the viscoelastic response alone, without the sea level equation |
| `kernel_int` | Integrates a kernel file over the model volume |
| `visc_corr` | Correlation and rms, as a function of radius, between two 3D viscosity fields |
| `merge_files`, `read_kl`, `read_kl_2`, `read_rob` | Post-processing of output files for plotting |


### Configuring a run

Every program has defaults built in, and takes an optional namelist file that
overrides them. So the simplest run is still

```
cd build/run
./sl3d/sea_level_3d prem.simp.vm1
```

and a configured one is

```
./sl3d/sea_level_3d run.nml
```

with `run.nml` holding as much or as little as you want to change:

```
&run
  earth_model    = 'prem.simp.vm1'
  ice_model_name = 'ice6g'      ! ice5g, ice6g or ice7g
  lmax_run       = 64           ! spherical harmonic truncation degree
  time_step      = 5.0          ! years
  output_every   = 50.0         ! years
  run_length     = -1.0         ! years; negative is the whole ice history
/

&kernel                         ! read by the kernel and inversion programs
  site_lat       = 57.0
  site_long      = -110.0
  load_degree    = 2
  load_order     = 2
/
```

Anything you leave out keeps the program's own default, so the same file can be
used with programs that want quite different settings. Both groups are
optional, as is every entry in them. Arguments are recognised by what they
name rather than by their position, so a config file, an earth model and an ice
model can be given in any combination and any order, and anything named on the
command line overrides the file:

```
./sl3d/sea_level_3d run.nml ice7g       # the file, but with a different ice model
./sl3d/sea_level_3d prem.simp.vm1 ice6g # no file at all
```

Each run writes `config.out.nml` beside its results, holding the parameters it
actually used. That file is itself a valid config file, so a run can be
repeated exactly by handing its record back to the program.


## Plotting

`scripts` holds a few small Python scripts for looking at forward modelling
output. They need only numpy and matplotlib, which most systems already have;
otherwise, from the top of the repository:

```
python3 -m venv venv && . venv/bin/activate
pip install -r scripts/requirements.txt
```

Nothing heavier is needed because the codes write plain text, and the map
projection is only twenty lines: `scripts/sl3dplot.py` reproduces the Robinson
projection from `module_gia.f90`, and the coastlines in
`data/coast_lines_110_rob.dat` were projected with it, so they land on the data
without any mapping library.

| Script | What it does |
| --- | --- |
| `plot_field.py` | A global map of one output file |
| `plot_rsl.py` | A sea level curve at one or more sites, from a run's output files |
| `plot_model.py` | Density, elastic moduli and viscosity against radius, from a model file |

Run them where the output is, which is normally `build/run`:

```
cd build/run
../../scripts/plot_field.py sl.1d.18.00
../../scripts/plot_field.py sl.1d.18.00 --column 3 -o sea_level.png
../../scripts/plot_rsl.py --lat 57.7 --lon 11.9 sl.1d.*
../../scripts/plot_model.py prem.simp.vm1
```

Each opens a window by default, or writes a png with `-o`, and `--help`
explains the rest.

The two codes write slightly different files: `gia1d` writes longitude,
latitude, sea level and its change since the start of the calculation, while
`sl3d` writes the output index and the projected coordinates as well. The
scripts work out which is which from the file, so the same commands serve both.
Columns are numbered as they are in the file, as they would be in gnuplot, and
`--column` defaults to the last, which in both cases is the change in sea
level; the column before it is the sea level itself.

`plot_rsl.py` takes the age of each output file from the end of its name, as
the forward programs write it, so `sl.1d.20.50` is read as 20.5 ka.
Give it a run's files and a site, and it assembles the sea level history there.
Sites can be named and repeated:

```
../../scripts/plot_rsl.py sl.1d.* \
    --site 57.7 11.9 Fennoscandia --site -17.5 -149.5 Tahiti
```

The kernel programs write their output through the spectral element mesh rather
than on a surface grid, so these scripts do not read it.

## Earth model files

A model file is plain text:

```
<name>
<nknot> <ntau>
r  rho  kappa  mu  log10(eta)
...                                (nknot rows)
```

with radius in m, density in kg m<sup>-3</sup>, bulk and shear moduli in GPa,
and the base 10 logarithm of the viscosity in Pa s. A viscosity entry of `0`
marks a region that is not viscous, i.e. the elastic lithosphere or the fluid
outer core. `ntau` is the number of Maxwell elements in the rheology, and is `1`
for a standard Maxwell solid. Knots are ordered outwards from the centre, and a
repeated radius marks a discontinuity between layers.

`data/prem.simp.vm1` is a simplified PREM elastic structure with a two layer
viscosity profile, 5x10<sup>20</sup> Pa s in the upper mantle over 5x10<sup>21</sup>
Pa s in the lower mantle, and can be used as a template.


## Other inputs

Besides the earth model given on the command line, the codes take these
inputs. The last two are still chosen in the source rather than at run time:

- **Ice history.** Three ice models are available, each held in `data` as one
  netCDF file per time, reduced from the standard distribution to the two
  fields the codes read plus the coordinates:

  | Name | Model | Span | Times |
  | --- | --- | --- | --- |
  | `ice5g` | ICE-5G v1.2 | 21 ka to present | 39 |
  | `ice6g` | ICE-6G_C VM5a | 26 ka to present | 48 |
  | `ice7g` | ICE-7G_NA VM7 | 26 ka to present | 48 |

  Choose one with an optional second argument, or with the `ICE_MODEL`
  environment variable; the default is `ice5g`:

  ```
  cd build/run
  ./sl3d/sea_level_3d prem.simp.vm1 ice6g
  ```

  Programs that integrate the whole history take their time window from the
  model, so selecting `ice6g` or `ice7g` also extends the calculation from
  21 ka to 26 ka at an unchanged step length. The files are read from a
  directory named after the model, which the build links into `build/run`;
  setting `ICE_DIR` looks in `$ICE_DIR/ice5g` and so on instead.

  The files are netCDF-3 classic and uncompressed, which git compresses to
  about half when packing, so compressing them here would save little.

  `module_ice.f90` holds everything that differs between the three, which is
  more than it looks: ICE5G names the ice thickness `sftgit` and ICE6G and
  ICE7G call it `stgit`; ICE5G's `orog` is the field that ICE6G and ICE7G call
  `Topo`, while ICE6G's own `orog` is clipped at zero and so carries no
  bathymetry; and ICE5G's longitudes sit on whole degrees where the other two
  sit on half degrees, so those two are averaged onto the whole degree grid
  that the codes interpolate on.
- **3D viscosity (sl3d).** `visc_3d_from_vs` in `module_model.f90` scales
  viscosity from the shear wave model `S20RTS.sph`, which is expected in the
  working directory; `visc_3d` sets perturbations directly in code; and
  `visc_3d_none` gives a laterally uniform model. Each program chooses one of
  these near the start of its source.
- **Resolution and time stepping.** The spherical harmonic truncation degree
  (`lmax`), the number of time steps (`nt`) and the step length (`dt`) are set
  as literals in each program, and changing them means recompiling.


## References
<a id="1">[1]</a>
Crawford, O., Al-Attar, D., Tromp, J., Mitrovica, J.X., Austermann, J. and Lau, H.C., 2018. Quantifying the sensitivity of post-glacial sea level change to laterally varying viscosity. Geophysical journal international, 214(2), pp.1324-1363.

<a id="2">[2]</a>
Crawford O., 2019. On the viscoelastic deformation of the Earth. Ph.D. thesis, University of Cambridge.
