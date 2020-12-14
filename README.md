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

- A function that takes in coordinates and gives you energies.

- A (minimized) geometry upon which to perform the harmonic analysis.

This code allows the user to perform 3 or 5 point finite difference for the first and second derivatives in the Hessian.  

For more information about finite difference, see:
http://www.mathematik.uni-dortmund.de/~kuzmin/cfdintro/lecture4.pdf
