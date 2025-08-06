function [t, x, hdata] = ercsj_scaled(dx, tspan, x0, D, ...
    Atol, Rtol, h, get_hdata, use_ercsj)
% Integrator parameters
buffer = 512;
threshold = Atol / Rtol;
t0 = tspan(1);
T = tspan(end);
hmin = 1e-12;
hmax = 0.1 * (T - t0);
rhmin = 0;
rhmax = 8;
p = 2; % Global order of estimator

% Initialize variables
n = length(x0);
I = eye(n, n);
if isempty(D)
    D = I;
    Di = I;
else
    Di = inv(D);
end

% Initial buffer of solution outputs
t = zeros(1, buffer);
x = zeros(n, buffer);
t(1) = t0;
x(:,1) = x0;

if isempty(h)
    % Estimate initial step
    % err1_h = max(abs(dx(t0, x0)) ./ max(abs(x0), threshold));
    % h = (0.8 ^ (p + 1)) * Rtol / err1_h;
    % h = max(hmin, min(hmax, h));
    h = initial_step(dx, x0, t0, Atol, Rtol, hmin, hmax);
    fixed_h = false;
else
    fixed_h = true;
end

if isempty(use_ercsj)
    use_ercsj = true;
end

hdata = [];
if get_hdata
% Step data is required
    hdata.hcurve = h;
    hdata.ht = t0;
    hdata.h_rej = false;
    hdata.e1 = zeros(n, 1);
    hdata.e2 = zeros(n, 1);
end

% Main loop
rejectStep = false;
i = 1;
done = false;
f_info = [];
while ~done
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
    
%     % Richardson extrapolation
%     if ~rejectStep || true
%         % One step of size h (Crude)
%         x_est = my_integrator_step(dx, x_i, t_i, h, D, Di, I, []);
%     end
%     % Two steps of size h/2 (Accurate)
%     x_mid = my_integrator_step(dx, x_i, t_i, h/2, D, Di, I, []);
%     x_new = my_integrator_step(dx, x_mid, t_i + h/2, h/2, D, Di, I, []);
%     % Error of the accurate step
%     e_i = (2^p) .* (x_new - x_est) ./ ((2^p) - 1);
%     x_new = x_est;
%     [~, e2_i] = my_integrator_step(dx, x_i, t_i, h, D, Di, I, []);

    % Integrate next step
    [x_new, e_i, f_info] = RCSJ_step(dx, x_i, t_i, h, D, Di, ...
        I, f_info, fixed_h);
    if use_ercsj
        x_new = x_new - e_i;
    end
    
    % Calculate relative error and accept/reject step
    err = max(abs(e_i) ./ max(max(abs(x_new), abs(x_i)), ...
        threshold));
    if (err > Rtol || isnan(err)) && ~fixed_h
        % Step rejected, halve step size and repeat step
        rejectStep = true;
        h = h / 2;
%         f_info = [];
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
        if t(i) >= T
            done = true; % Last step accepted
        end

        if ~fixed_h
            % Update step size
            rh = 0.8 / ((err / Rtol) ^ (1 / (p + 1)));
            rh = max(rhmin, min(rhmax, rh));
            h = max(hmin, min(hmax, h * rh));
        end
    end
    
    if get_hdata
        % Save new step size in curve.
        hdata.hcurve = [hdata.hcurve h];
        if rejectStep
            hdata.ht = [hdata.ht t_i];
            hdata.h_rej(end) = true;
        else
            hdata.ht = [hdata.ht t(i)];
        end
        hdata.h_rej = [hdata.h_rej false];
        hdata.e1 = [hdata.e1 e_i];
        if exist('e2_i', 'var')
            hdata.e2 = [hdata.e2 -e2_i];
        else
            hdata.e2 = [hdata.e2 e_i];
        end
    end
    
    % If the last step was rejected and the step size is too small, the
    % integration failed.
    if rejectStep && h < hmin
        warning(['The integration failed to achieve the specified ' ...
            'tolerance at time ' num2str(t_i) '. The solution ' ...
            'up to this time will be returned.']);
        break;
    end
end
% Trim output vectors
t = t(1:i);
x = x(:,1:i);


function [x_new, err, f_info] = RCSJ_step(f, x0, t0, h, D, ...
    Di, I, f_info, fixed_h)
% A-stable, order 2 step
c = [1 -4/3 -5/3 2/3];
cGi = [-7/6 3 1/3 -5/6];
n = length(x0);
F = zeros(n, n);
for j = 1:n
    F(:,j) = f(t0, x0 + 1i .* 0.5 .* h .* D(:,j));
end
FI = imag(F * Di);
[L,U,p] = lu(I - FI, 'vector');
if ~fixed_h
    adU = abs(diag(U));
    tol = max(size(U)) * eps(max(adU)); % taken from rank.m
    if any(adU <= tol)
        % Singular matrix, reject step
        x_new = nan(size(x0));
        err = nan(size(x0));
        f_info = [];
        return;
    end
end

if isempty(f_info)
    v0 = f(t0, x0);
    v3 = f(t0 + h, x0);
else
    v0 = f_info(:,1);
    v3 = f_info(:,2);
end
v1 = f(t0 + 0.5 .* h, x0);
Giv1 = U \ (L \ v1(p));
x_new = x0 + h .* Giv1;

if nargout > 1
    % Evaluate f at the required points
    v2 = f(t0 + 0.5 .* h, x0 + 0.5 .* h .* Giv1);
    
    % Asymptotic error term
    f_sum = cGi(1).*v0 + cGi(2).*v1 + cGi(3).*v2 + cGi(4).*v3;
    f_sum = U \ (L \ f_sum(p));
    f_sum = f_sum + c(1).*v0 + c(2).*v1 + c(3).*v2 + c(4).*v3;
    err = h .* f_sum;
    
    % Stabilization
    err = U \ (L \ err(p));
    err = err - 1.5 .* (FI*err);
    err = U \ (L \ err(p));
    
    % Update info of known f values
    f_info = [v0 v1];
end