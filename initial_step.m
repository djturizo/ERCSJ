function h = initial_step(dx, x0, t0, Atol, Rtol, hmin, hmax)
% Initial step size is calculated as the optimal step for the Euler method 
% at the initial point, using Richardson extrapolation for error estimation

threshold = Atol / Rtol;
p = 1; % Euler method is of order 1

dx0 = dx(t0, x0);
% Calculate small step for extrapolation
h0 = 0.01 * max(abs(x0)) / max(abs(dx0) ./ max(abs(x0), threshold));
h0 = max(hmin, min(hmax, h0));
x1_est = x0 + h0 .* dx0;
x1 = x0 + (h0/2) .* dx0;
x1 = x1 + (h0/2) .* dx(t0 + h0 / 2, x1);
% e1 = (x1 - x1_est) / ((2^p) - 1);
e1 = (x1_est - x1) / ((2^p) - 1);
err = max(abs(e1) ./ max(max(abs(x1), abs(x0)), threshold));
% rh can tend to inf if err is too small, better use irh = 1/rh
irh = ((err / Rtol) ^ (1 / (p + 1))) / 0.8;
if h0 > hmax * irh;
    h = hmax;
else
    h = max(hmin, min(hmax, h0 / irh));
end

