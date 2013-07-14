math-tools
==========


###Complex Functions###
This folder contains three scripts, complexanalysis.py, complexfunctions.py, and complexgrapher.py.
I wrote these scripts after being inspired by a Complex Analysis course taught at Andover (Math 630).

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
computes the images of horizontal and vertical lines under the given function. These
graphs are based off the ideas mentioned in Visual Complex Analysis by Tristan Needham.

Examples of graphs:

#####Sine#####
![](complex-functions/sine.png?raw=true)

Note the real period of 2pi.

#####Exponential#####
![](complex-functions/exp.png?raw=true)

Note the imaginary period of 2pi*i

#####Logarithm#####
![](complex-functions/log.png?raw=true)

Note the branch cut along the real axis.

#####Gamma#####
![](complex-functions/gamma.png?raw=true)

Note the singularities at the non-positive integers. These singularities are of order 1,
since the colors wrap around these singularities only once.

#####Elliptic Theta#####
![](complex-functions/Theta.png?raw=true)

Note the real period. Also, although the theta function is not periodic along the imaginary axis,
all the zeros line up in the imaginary direction.

#####sn#####
![](complex-functions/sn.png?raw=true)

Note the period along the real axis and that along the imaginary axis.
Also, notice that each "cell" has a singularity and a zero -- as predicted by Liouville's theorem.

#####log square plot#####
![](complex-functions/log_square.png?raw=true)

Notice all the sheets, corresponding to the Riemann surface of the log function.
This is once again related to the branch cut along the negative real axis.

#####sn square plot#####
![](complex-functions/sn_square.png?raw=true)

Notice once again the double-periods: the real lines and imaginary lines get mapped to 
closed loops -- except those lines that get mapped to infinity.
