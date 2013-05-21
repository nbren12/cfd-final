# README
## Summary

I implemented the numerical schemes using a combination of FORTRAN and python. FORTRAN was used for the computationally expensive reimann solver, boundary condition treatment, WENO reconstruction. Python was used "glue". A Controller class was implemented that among other things, writes the data to disk, loops through the time stepping, and passes the required parameters to the FORTRAN scripts.

The most important files are contained in the ./src/swe2d directory. `hr_solver2d.f90` and `rp_roe.f90` contain the High Resolution solver and roe averaged riemann solver.

The code for the radial dam break, 1d dam break, and geostrophic adjustment problems is contained in the examples directory. Executing `python examples/rad-dam-break/rad-dam-2d.py` should make a directory `n100_hr` containing saved output. This folder can be imported into a python shell by `import n100_hr`. Then, `n100_hr.cont` is a instance of the ControllerSW2D class. `cont.state` is an instance of the `State` class defined in `pyclaw`. `n100_hr.cont.state.q` contains the data.

Compilation and packaging is by accomplished `numpy.distutils`.


## Directory Structure

    ├── README.md            -- this README file
    ├── examples             -- folder with examples
    │   ├── adv-weno.py      -- simple weno advection example
    │   ├── dam-break        -- contains code for the 1d dam break 
    │   ├── geo-adjustment   -- contains code for the geostrophic adjustment
    │   └── rad-dam-break    -- contains code for the radial dam break
    ├── src
    ├── swe1d                -- contains old 1D implementation
    ├── swe2d                -- contains the 2D SWE implementation
    │   ├── py_solvers.py    -- contains the source terms solver for the radial dam break
    │   ├── hr_solver2d.f90  -- time stepper that drives the rp_roe.f90 riemann solver
    │   |                       Also contains the time stepper for the coriolis terms.
    │   ├── Controller.py    -- the custom Controller class
    │   ├── bc2d.f90         -- boundary conditions
    │   ├── rp_roe.f90       -- Riemann solver,flux corrector, and limiters
    └── weno                 -- contains code for WENO schemes
        ├── adv_weno.f90     -- the time advancing 
        └── weno.f90         -- The WENO reconstruction. Courtesy of PyWENO.

## Dependencies

- numpy
- scipy
- matplotlib
- [clawpack](https://github.com/clawpack/clawpack). Just used for some convenience classes. I implemeneted the meat of the code myself.

## Installation

Install the dependencies. `clawpack` can be installed with `pip install clawpack`. For the other dependecies, I recommend using some standard python distribution.

Checkout the code: `git clone https://github.com/nbren12/cfd-final.git cfd-final`

Change to `cfd-final` and either do (recommended):


        python setup.py build_ext --inplace
        export PYTHONPATH=<path to cfd-final>/src:$PYTHONPATH

or (not recommended, but probably works)

        python setup.py install

Then go to the `examples` directory and poke around. The code should work if the installation was succesful.
