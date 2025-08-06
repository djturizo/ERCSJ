function [t, x] = ercsj(f, tspan, x0, varargin)
%'ercsj' Solve stiff differential equations that admit a holomorphic 
% extension.
%   [t_sol, x_sol] = ercsj(f, tspan, x0)
%   Solves the differential equation x' = f(t,x). 'f' is a function
%   handle with 2 arguments: 't' as the time and 'x' as the state variable,
%   in that order. 'f' must return a column vector or a scalar. 
%   'tspan'=[t0 Tf] defines the integration interval. The initial condition
%   is given by x(t0) = x0. The numerical solution of the state variable at
%   the times given by 't_sol' is 'x_sol'.
%
%   [t_sol, x_sol] = ercsj(f, tspan, x0, Rtol)
%   [t_sol, x_sol] = ercsj(f, tspan, x0, Rtol, Atol)
%   Rtol is the relative error tolerance (1e-3 by default) of the numerical
%   method. Atol is the absolute error tolerance (1e-6 by default).

% Default parameters
narginchk(3,5);
Rtol = 1e-3;
Atol = 1e-6;
if nargin == 4
    Rtol = varargin{1};
elseif nargin == 5
    Atol = varargin{2};
end

% Integrator parameters
buffer = 512;
threshold = Atol / Rtol;
t0 = tspan(1);
T = tspan(end);
hmin = 1e-12;
hmax = 0.1 * (T - t0);
rhmin = 0;
rhmax = 8;
order = 2; % Global order of estimator

% Equation constants
c = [1 -4/3 -5/3 2/3];
cGi = [-7/6 3 1/3 -5/6];

% Initialize variables
n = length(x0);
I = eye(n, n);
FI = zeros(n, n);

% Initial buffer of solution outputs
t = zeros(1, buffer);
x = zeros(n, buffer);
t(1) = t0;
x(:,1) = x0;

% Estimate initial step
h = initial_step(f, x0, t0, Atol, Rtol, hmin, hmax);

% Main loop
rejectStep = false;
i = 1;
f_info = [];
while true
    % If the last step was rejected and the step size is too small, the
    % integration failed.
    if rejectStep && h < hmin
        warning(['The integration failed to achieve the specified ' ...
            'tolerance at time ' num2str(t_i) '. The solution ' ...
            'up to this time will be returned.']);
        break;
    end
    
    if ~rejectStep
        % Advance timestep
        t_i = t(i);
        x_i = x(:,i);
        i = i + 1;
    end
    
    % Avoid exiting time span.
    if t_i + h > T
        h = T - t_i;
    end
    
    % Integrate next step (A-stable, order 2 step)
    for j = 1:n
        FI(:,j) = imag(f(t_i, x_i + 1i .* 0.5 .* h .* I(:,j)));
    end

    [L,U,p] = lu(I - FI, 'vector');
    adU = abs(diag(U));
    tol = max(size(U)) * eps(max(adU)); % taken from rank.m
    if any(adU <= tol)
        % Singular matrix, reject step and halve step size
        rejectStep = true;
        h = h / 2;
        f_info = [];
        continue;
    else
        % Evaluate f at the required points
        if isempty(f_info)
            v0 = f(t_i, x_i);
            v3 = f(t_i + h, x_i);
        else
            v0 = f_info(:,1);
            v3 = f_info(:,2);
        end
        v1 = f(t_i + 0.5 .* h, x_i);
        Giv1 = U \ (L \ v1(p));
        x_new = x_i + h .* Giv1;

        v2 = f(t_i + 0.5 .* h, x_i + 0.5 .* h .* Giv1);

        % Asymptotic error term
        f_sum = cGi(1).*v0 + cGi(2).*v1 + cGi(3).*v2 + cGi(4).*v3;
        f_sum = U \ (L \ f_sum(p));
        f_sum = f_sum + c(1).*v0 + c(2).*v1 + c(3).*v2 + c(4).*v3;
        e_i = h .* f_sum;
        
        % Stabilization
        e_i = U \ (L \ e_i(p));
        e_i = e_i - 1.5 .* (FI*e_i);
        e_i = U \ (L \ e_i(p));

        % Update info of known f values
        f_info = [v0 v1];
    end
    
    % ELICS correction
    x_new = x_new - e_i;
    
    % Calculate relative error and accept/reject step
    err = max(abs(e_i) ./ max(max(abs(x_new), abs(x_i)), ...
        threshold));
    if err > Rtol || isnan(err)
        % Step rejected, halve step size and repeat step
        rejectStep = true;
        h = h / 2;
    else
        % Step accepted
        rejectStep = false;
        f_info = [];

        % Save step
        if i > length(t)
            % Allocate more memory to output vectors
            t = [t zeros(1, buffer)];
            x = [x zeros(n, buffer)];
        end
        t(i) = t_i + h;
        x(:,i) = x_new;
        
        % Break loop when the integrstion finishes
        if t(i) >= T
            break;
        end

        % Update step size
        rh = 0.8 / ((err / Rtol) ^ (1 / (order + 1)));
        rh = max(rhmin, min(rhmax, rh));
        h = max(hmin, min(hmax, h * rh));
    end
end
% Trim output vectors
t = t(1:i);
x = x(:,1:i);
