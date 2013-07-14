math-tools
==========


###Complex Functions###
This folder contains three scripts, complexanalysis.py, complexfunctions.py, and complexgrapher.py

complexanalysis.py has some code for computing path integrals, residues, Laurent coefficients, etc.
It is meant more as a demonstration rather than an actual computational library.

complexfunctions.py has some code defining several complex functions, like the gamma function, the 
Riemann zeta function, the Weierstrass Elliptic function, and more. Once again, no claims are made about
computational efficiency; it is for demonstrations too.

complexgrapher.py has two main functions for plotting complex functions, `plot_function` and `square_plot`.
`plot_function` computes the values of functions on a specified square grid, and colors the points on the
grid according to the function output: the modulus is indicated by shade (larger modulus correspond to
brighter colors), and the argument is indicated by color (the real axis is red, numbers with argument
of 120 degrees are green, and those with 240 degrees are blue). `square_plot`, on the other hand,
computes the images of horizontal and vertical lines under the given function.
