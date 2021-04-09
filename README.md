# Optimization-Based Relaxations of Parametric ODEs
Consider the following parametric system of ordinary differential equations (ODEs). Here `t` denotes time, `p` denotes a vector of system parameters,  `x` denotes a vector of system states, and `x'` is the time-derivative of `x`.
```
x' = f(t,p,x),     x(0,p) = x0(p)
```
This repository contains a proof-by-concept implementation in MATLAB of a [new method by Song and Khan]() to construct and evaluate useful state relaxations for this ODE system. Roughly, given `f` and `x0`, a final time `tf`, and parameter bounds `[pL,pU]`, this implementation constructs convex and concave relaxations of each component of `x(tf,p)` with respect to `p` on `[pL,pU]`, by generating and solving an auxiliary ODE system. These relaxations enclose the reachable set of the original ODE system, and may be used in deterministic methods for continuous global dynamic optimization.

This implementation was developed by Yingkai Song.

# Method outline
to be written

# Usage
to be written

# References
to be written
