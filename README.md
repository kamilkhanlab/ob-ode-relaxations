# Optimization-Based Relaxations of Parametric ODEs
Consider the following parametric system of ordinary differential equations (ODEs). Here *t* denotes time, **p** denotes a vector of system parameters, **x** denotes a vector of system states, and **x**' is the time-derivative of **x**.

> **x**'(*t*,**p**) = **f**(*t*,**p**,**x**(*t*,**p**)),
> 
> **x**(0,**p**) = **x**<sub>0</sub>(**p**)

This repository contains a proof-by-concept implementation in MATLAB of a [new method by Song and Khan]() to construct and evaluate useful state relaxations for this ODE system. Roughly, given functions **f** and **x**<sub>0</sub>, a final time *t*<sub>f</sub>, and parameter bounds [**p**<sup>L</sup>, **p**<sup>U</sup>], this implementation constructs convex underestimators and concave overestimators of each component of the final state **x**(*t*<sub>f</sub>,**p**) with respect to **p** on [**p**<sup>L</sup>, **p**<sup>U</sup>], by generating and solving an auxiliary ODE system. These relaxations enclose the reachable set of the original ODE system, and may be used in deterministic methods for continuous global dynamic optimization.

This implementation was developed by Yingkai Song.

# Method outline
to be written

# Usage
to be written

# References

- Y. Song and K.A. Khan, Optimization-based convex relaxations for nonconvex parametric systems of ordinary differential equations, *Math Program*, accepted.
