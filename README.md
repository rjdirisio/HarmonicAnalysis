# HarmonicAnalysis
Python code to do harmonic analysis of arbitrary molecules

This Harmonic Analysis code generates a Hessian Matrix, Mass-weights it, and then diagonalizes it to obtain squares of
harmonic frequencies and mass-weighted eigenvectors (vibrational normal modes).

What this code requires:
-A function that takes in coordinates and gives you energies.
-A geometry upon which to perform the harmonic analysis.

This code, at the moment, performs a 5-point finite difference for the second derivatives in the Hessian.  It also uses 4 points out of a 9 point stencil for the mixed derivatives.  For more information about this, see:
http://www.mathematik.uni-dortmund.de/~kuzmin/cfdintro/lecture4.pdf

I will perhaps add in flexibility of higher order differences later on.
