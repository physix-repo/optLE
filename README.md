# Langevin Equation Optimization

This program optimizes the parameters of a Langevin equation to fit trajectories
from molecular dynamics.

## Prerequisites

* FORTRAN compiler

## Installing

Run `make` inside the project root.

The objects will be generated in `build`, and the executables will be generated
in `bin` directory.

## Usage

This program reads three input files: _input_, _colvar_ and _RESTART_. All these
files have to be in the same directory. A detailed description of the files is
shown below.

After compilation, you can run the program inside the `bin` directory:

```
cd bin && ./optle
```

### _input_ file

This is where the program parameters are stored in.

An standard input file is copied into the bin dir during the `make` task. This
file contains all the parameters and their description. This file is also
available in [resources/input](./resources/input).

The input file is the following:

```
colvar_file       colvar        (format: t x, multiple trajectoires appended)  
dt                0.01          (time step in colvar file)
dtmult            1             (time step for integration dtint = dt/dtmult)
xmin            -3.0            (min and max value of the CV in the colvar)
xmax             3.0
xbins            40             (used for x and v histograms in old version: not used now)
kT               1.0
init_tau         400.           (not used in overdamped optimization)
opt_niter        100000         (optimization steps, 0 = just print one traj and stop)
opt_temp         3e-4 3e-4 0.05    (annealing T: initial, final, target acceptance)
type_Langevin      0            (0 = overdamped, 1 = std, 2 = GLE)
ratio_Langevin_MD  1            (ratio between number of Langevin traj and MD traj)
fit_F              1            (free energy profile fit 1 = on, 0 = off)
fit_gamma          1            (friction profile fit 1 = on, 0 = off)
fit_tau            0
fit_mass           0            (mass profile fit 1 = on, 0 = off)
type_error         3            (1 = -logLik(KL), 2 = RMSD, 3 = -logLik(analyt.propagator (-3=test&exit)), 4 = -logLik(num.propagator))
pos_dep_gamma      1            (0 = fixed gamma, 1 = position-dependent gamma)
pos_dep_mass       0            (0 = fixed mass , 1 = position-dependent mass )
max_Gaussian_h   2.  2.   1.   (max height of Gaussians added to profiles F,g,m)
max_Gaussian_w   0.05 0.05 0.05 (max width of Gaussians added, units of (xmax-xmin))
fix_mass0          0            (keep fixed the value of the mass at the TS)
use_velocity       0            (0 = optimize P(x,t), 1 = optimize P(x,v,t))
dtint_prop         0.1          (time step for integration in propagator test)
ntraj_prop         10           (number of trajs to test the propagator) 
```

### _colvar_ file

This file contains multiple trajectories appended, in the following format:

`time (t)`, `collective variable (x)`

Example:

```
0.00000000E+00 0.00000E+00
0.10000000E-01 -0.15808E-01
0.20000000E-01 -0.14025E-01
0.30000000E-01 -0.79028E-05
0.40000000E-01 0.16491E-01
…
0.99800000E+01 0.10768E+01
0.99900000E+01 0.10805E+01
0.10000000E+02 0.10787E+01
```

### _RESTART_ file

This file contains an initial guess for the free energy (F), friction (gamma)
and mass profiles. It must be five columns in this file, separated by spaces:

```
# x F F/kT gamma mass
0.0580   0.0000   0.0000 200.0000   1.0000
0.0600   0.0000   0.0000 200.0000   1.0000
0.0620   0.0000   0.0000 200.0000   1.0000
0.0640   0.0000   0.0000 200.0000   1.0000
0.0660   0.0000   0.0000 200.0000   1.0000
0.0680   0.0000   0.0000 200.0000   1.0000
…
```

### Output

The important output files of the program are:

|File     |Description                                          |
|---------|-----------------------------------------------------|
|PROFILES |F(q) and D(q) for the last accepted optimization step|
|fort.5***|F(q) and D(q) of any accepted optimization step      |
|log      |Optimization status, -log(L), initial D estimate, ...|
