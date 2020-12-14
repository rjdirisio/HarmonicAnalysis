# HarmonicAnalysis
Minimal Python package to do harmonic vibrational normal mode analysis of arbitrary molecules

## Installation

To install the package, `git clone` the directory from GitHub, then `cd` into the directory.
After, you may run

`pip install .`

In order to install the package, which can then be used in an `import`

## Description

This Harmonic Analysis code generates a Hessian Matrix, Mass-weights it, and then diagonalizes it to obtain 
harmonic frequencies and mass-weighted eigenvectors (vibrational normal modes).

What this code requires:

- PyVibDMC. This code uses the potential_manager outlined in the PyVibDMC documentation

- A function that takes in coordinates and gives you energies (set up with potential_manager)

- A (minimized) geometry upon which to perform the harmonic analysis.

This code allows the user to perform 3 or 5 point finite difference for the first and second derivatives in the Hessian.  

For more information about finite difference, see:
http://www.mathematik.uni-dortmund.de/~kuzmin/cfdintro/lecture4.pdf


## Example

```
from HarmonicAnalysis import * # Imports harmonic_analysis

# Everything is in  Atomic Units going into generating the Hessian.
dxx = 1.e-3
water_geom = np.array([[0.9578400, 0.0000000, 0.0000000],
                       [-0.2399535, 0.9272970, 0.0000000],
                       [0.0000000, 0.0000000, 0.0000000]])

#Specify potential manager
pot_dir = 'path/to/Partridge_Schwenke_H2O/'
py_file = 'h2o_potential.py'
pot_func = 'water_pot'
partridge_schwenke = pm.Potential(potential_function=pot_func,
                                  potential_directory=pot_dir,
                                  python_file=py_file,
                                  num_cores=1)

# Using pyvibdmc.simulation_utilities
geom = Constants.convert(water_geom, "angstroms", to_AU=True)  # To Bohr from angstroms

#Specify atoms
atms = ["H", "H", "O"]

harm_h2o = harmonic_analysis(eq_geom=geom,
                          atoms=atms,
                          potential=partridge_schwenke,
                          dx=dxx)
freqs, normal_modes = harm_h2o.run()
# Turns off scientific notation
np.set_printoptions(suppress=True)
print(f"Freqs (cm-1): {freqs}")
```
