# ERCSJ
Extrapolated Rosebrock method with complex step Jacobian estimator for numerical Integration of stiff ordinary differential equation that admit a holomorphic extension.

## Syntax

    [t_sol, x_sol] = ercsj(f, tspan, x0);
    
Solves the differential equation $x' = f(t,x)$. `f` is a function handle with 2 arguments: $t$ as the time and $x$ as the state variable, in that order. $f$ must return a column vector or a scalar.  `tspan`$=[t_0, T_f]$ defines the integration interval. The initial condition is given by $x(t_0) =$ `x0`. The numerical solution of the state variable at the times given by `t_sol` is `x_sol`.

    [t_sol, x_sol] = ercsj(f, tspan, x0, Rtol);
    [t_sol, x_sol] = ercsj(f, tspan, x0, Rtol, Atol);
    
`Rtol` is the relative error tolerance (`1e-3` by default) of the numerical method. `Atol` is the absolute error tolerance (`1e-6` by default).
