# ob-ode-relaxations

Consider the following parametric system of ordinary differential equations (ODEs). Here $t$ denotes time, $\mathbf{p}$ denotes a vector of system parameters, $\mathbf{x}$ denotes a vector of system states, and $\dot{\mathbf{x}}$ is the time-derivative of $\mathbf{x}$.

$$
\begin{align*}
\dot{\mathbf{x}}(t,\mathbf{p}) &= \mathbf{f}(t,\mathbf{p},\mathbf{x}(t,\mathbf{p})) \\
\mathbf{x}(t_0,\mathbf{p}) &= \mathbf{x}_0(\mathbf{p}).
\end{align*}
$$

This repository contains a proof-by-concept implementation in MATLAB of a [new method by Song and Khan][1] to construct and evaluate useful state relaxations for the ODE system above. Roughly, given functions **f** and **x**<sub>0</sub>, a time horizon [*t*<sub>0</sub>,*t*<sub>f</sub>], and parameter bounds [**p**<sup>L</sup>, **p**<sup>U</sup>], this implementation constructs convex underestimators and concave overestimators of each component of the final state **x**(*t*<sub>f</sub>,**p**) with respect to **p** on [**p**<sup>L</sup>, **p**<sup>U</sup>], by generating and solving an auxiliary ODE system. These relaxations enclose the reachable set of the original ODE system, and may be used in deterministic methods for continuous global dynamic optimization.

The following components of this implementation may be used on their own:

- a straightforward MATLAB implementation of interval arithmetic using operator overloading. This automatically computes upper and lower bounds for a composite function (such as **f** above) on a box-shaped domain.
- a straightforward MATLAB implementation of generalized McCormick relaxation using operator overloading, analogous to [MC++](https://github.com/coin-or/MCpp). This automatically computes convex and concave relaxations of each component of a composite function (such as **f** above) on a box-shaped domain.

This implementation was developed by Yingkai Song, and was used in all of the numerical examples in the [accompanying paper][1] except Example 6. This repository is tied to the accompanying paper, and will not be updated except for bug fixes. If you use this implementation, please cite our article:

> Y. Song and K.A. Khan, Optimization-based convex relaxations for nonconvex parametric systems of ordinary differential equations, *Math Program*, **196**, 521-565 (2022). doi:10.1007/s10107-021-01654-x

This work was supported by the Natural Sciences and Engineering Research Council of Canada (NSERC) under Grant RGPIN-2017-05944.

## Method outline
Interval arithmetic proceeds by replacing each arithmetic operation with corresponding bounding rules (see Moore (1979)), and generalized McCormick relaxation (Scott et al. (2011)) proceeds by appending additional convex relaxation rules to these bounding rules.

A new method by Song and Khan (2021) computes novel ODE relaxations by constructing and solving an auxiliary parametric ODE system with embedded convex optimization problems, whose objective functions employ convex and concave relaxations of the original right-hand side **f**. These optimal-value functions replace the flattened generalized McCormick relaxations **u**/**o** used in the Scott-Barton ODE relaxations (2013).

These new ODE relaxations require **p**–independent state bounds for **x**, which are constructed automatically via operator overloading using Harrison's bounding method (1977). Corresponding relaxations of **f** are constructed automatically via operator overloading using the McCormick relaxation method (McCormick (1976)). The new ODE relaxations are guaranteed to be at least as tight as the Scott–Barton relaxations (2013). Refer to [Song and Khan][1] for more details.  

## Usage

The files [example1.m](test/example1.m) and [example2.m](test/example2.m) in the [test](test) folder demonstrate the usage of this implementation. This implementation was developed and tested using MATLAB v9.6.0.1072779 (R2019a). It requires the `fmincon` solver in MATLAB's Optimization Toolbox. Each evaluation of the new relaxation ODE right-hand sides **u**/**o** requires solving convex optimization problems with `fmincon`; a more careful implementation would avoid this by adopting an active-set approach instead.

### Required inputs

This implementation applies to any parametric ODEs with factorable initial-condition and right-hand side functions. Supported operators in these functions are `-`, `+`, `.*`, `./`, `exp`, `log`, and `x.^n` where `n` is a positive integer between 2 and 21. This implementation also allows user-chosen MATLAB ODE solvers and tolerance settings. The following quantities and functions are required from users:

* `pL,pU`: lower and upper bounds of parameters **p**
* `p`: parameter vector at which relaxations will be computed
* `tspan`: time horizon `[t0,tf]` of the original parametric ODEs
* `original_initial_value(p)`: initial-condition function **x**<sub>0</sub> in the original parametric ODE
* `original_RHS(t,p,x,i)`: i<sup>th</sup> component of the original right-hand side function **f**
* `ODE_solver`: user-chosen MATLAB ODE solver
* `ODE_solver_options`: user-defined MATLAB ODE solver options
* `fmincon_options`: user-defined MATLAB options for the NLP solver `fmincon`

#### Caution

* Defining the functions `original_initial_value` and `original_RHS` above should only involve element-wise operators, such as `.*`, `./`, and `.^`. Matrix operators such as `*`, `/`, and `^` are not supported.
* If `original_initial_value` is independent of `p`, e.g., `original_initial_value(p)=0.82`, then trivial p-dependence must be added, e.g., `original_initial_value(p)=0.82+0.*p`. This allows MATLAB's operator overloading to access the initial condition.
* This implementation relies on valid Harrison state bounds. Thus, if Harrison bounds explode for the given parametric ODEs, this implementation cannot yield valid state relaxations. 

### Computing state relaxations

Add all the contents in the [src](src) folder to the MATLAB current folder. Then,
given all the required inputs above, run 

     [t,xAug] = compute_state_relaxations(p,pL,pU,tspan,@original_initial_value,@original_RHS,ODE_solver,ODE_solver_options,fmincon_options)
Outputs:

* `t`: a vector of time steps employed by the ODE solver for solving the auxiliary ODE system
* `xAug`: a matrix containing Harrison state bounds **x**<sup>L</sup>/**x**<sup>U</sup> of **x** on `[pL,pU]` at each time `t(i)`, and new state relaxations **x**<sup>cv</sup>/**x**<sup>cc</sup> of **x** on `[pL,pU]` at parameter `p` at each time `t(i)`. Row `i` of `xAug` contains the following quantities at time `t(i)` and parameter `p`:

$$
[x^{\mathrm{L}}_1, \ldots, x^{\mathrm{L}}_n, 
x^{\mathrm{U}}_1, \ldots, x^{\mathrm{U}}_n,
x^{\mathrm{cv}}_1, \ldots, x^{\mathrm{cv}}_n,
x^{\mathrm{cc}}_1, \ldots, x^{\mathrm{cc}}_n].
$$

## Implementation contents

This section outlines the various MATLAB scripts in the [src](src) folder.

### [Interval.m](src/Interval.m)

This class provides an implementation of natural interval extension (Moore (1979)) for factorable functions via operator overloading. Currently supported operators: `-`, `+`, `.*`, `./`, `exp`, `log`, and `x.^n` where `n` is a positive integer. 

This implementation can be used on its own. For example, if we have a MATLAB function of the form:
```Matlab
function y = f(x)
    % compute y using some combination of the supported operations
return
```
where `x` is a vector and `y` is a scalar, then the following code computes the natural interval extension `yI:=[yL,yU]` of `f` at an interval `xI:=[xL,xU]`. 
```Matlab
% express xL and xU as a single vector xI of Interval objects
nX = length(xL);
xI(1:nX) = Interval(0,0);
for i = 1:nX
    xI(i) = Interval(xL(i), xU(i));
end

% evaluate natural interval extension of f at xI
yI = f(xI);

% retrieve output bounds yL and yU, for which yL<=f(x)<=yU whenever xL<=x<=xU
yL = yI.lower;
yU = yI.upper;
```

### [McCormick.m](src/McCormick.m)

This class provides an implementation of generalized McCormick relaxations (Scott et al. (2011)) for factorable functions via operator overloading. Currently supported operators: `-`, `+`, `.*`, `./`, `exp`, `log`, and `x.^n` where `n` is a positive integer between 2 and 21.

This implementation can be used on its own. For example, if we have a MATLAB function of the form:
```Matlab
function y = f(x)
    % compute y using some combination of the supported operations
return
```
where `x` is a vector and `y` is a scalar, then the following code evaluates McCormick convex and concave relaxations of `f` at a point `x` in a box domain `[xL,xU]`. 
```Matlab
% express input data as a vector xMC of McCormick objects
nX = length(x);
xMC(1:nX) = McCormick(0,0,0,0);
for i = 1:nX
    xMC(i) = McCormick(xL(i), xU(i), x(i), x(i));
end

% evaluate McCormick relaxations of f at x on the domain [xL,xU]
yMC = f(xMC);

% retrieve corresponding convex and concave relaxations of f evaluated at x
yCv = yMC.convex;
yCc = yMC.concave;
```


### [convex\_relaxation\_of\_original\_RHS.m](src/convex_relaxation_of_original_RHS.m)  

A function `convex_relaxation_of_original_RHS(t,p,x,xL,xU,pL,pU,i,original_RHS)`, which evaluates McCormick convex relaxations of `original_RHS(t,p,x,i)` with `xL<=x<=xU` and `pL<=p<=pU` (see Assumption 3 in the article). 

### [concave\_relaxation\_of\_original\_RHS.m](src/concave_relaxation_of_original_RHS.m)  

A function `concave_relaxation_of_original_RHS(t,p,x,xL,xU,pL,pU,i,original_RHS)`, which evaluates McCormick concave relaxations of `original_RHS(t,p,x,i)` with `xL<=x<=xU` and `pL<=p<=pU` (see Assumption 3 in the article). Users are also allowed to replace these McCormick relaxations by other user-defined relaxations. 

### [linear\_transformation.m](src/linear_transformation.m)

A function `linear_transformation(alpha,xicv,xicc)`, which performs the linear transformation **v** described in equation (7) in the article. 

### [optimization\_based\_ODE\_RHS.m](src/optimization_based_ODE_RHS.m)

A function `optimization_based_ODE_RHS(t,xAug,p,pL,pU,original_RHS,fmincon_options)`, that constructs the RHS functions **u**/**o** of the new auxiliary ODE system (4) with embedded (8) in the article.

### [compute\_state_relaxations.m](src/compute_state_relaxations.m)

A function 

    compute_state_relaxations(p,pL,pU,tspan,original_initial_value,original_RHS,ODE_solver,ODE_solver_options,fmincon_options)
that solves the new auxiliary ODE system, yielding valid Harrison state bounds and new state relaxations.

## References

- Y. Song and K.A. Khan, [Optimization-based convex relaxations for nonconvex parametric systems of ordinary differential equations][1], *Math Program*, **196**, 521-565 (2022).
- J.K. Scott and P.I. Barton, Improved relaxations for the parametric solutions of ODEs using differential inequalities, *J Glob Optim*, **57**(1), 143-176 (2013)
- J.K. Scott, M.D. Stuber, and P.I. Barton, Generalized McCormick relaxations, *J Glob Optim*, **51**(4), 569-606 (2011)
- G. Harrison, Dynamic models with uncertain parameters. In: Avula, X. (eds.) Proceedings of the 1<sup>st</sup> International Conference on Mathematical Modeling, **1**, 295-304 (1977)
- R.E. Moore, *Methods and Applications of Interval Analysis*, SIAM, Philadelphia (1979)
- G.P. McCormick, Computability of global solutions to factorable nonconvex programs: Part I–Convex underestimating problems, *Math Program*, **10**(1), 147-175 (1976)

[1]: https://link.springer.com/article/10.1007/s10107-021-01654-x
