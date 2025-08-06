function test_ERCSJ(varargin)
close all;
case_num = 1;
fixed_h = false;
get_hdata = true;
do_comparison = false;
if nargin >= 1
    case_num = varargin{1};
end
if nargin >= 2
    fixed_h = varargin{2};
end
if nargin >= 3
    get_hdata = varargin{3};
end
if nargin >= 4
    do_comparison = varargin{4};
end
if fixed_h
    get_hdata = false;
end

if case_num == 1
    h = 1;
    t0 = 0;
    tf = 30;
    % A = [0 1; -c -b]; l^2+b*l+c=0;
    l1 = -1 + 1i.*0.8; l2 = conj(l1);
    % l1 = 1i.*0.8; l2 = conj(l1);
    % l1 = -1.5; l2 = -0.5;
    a = real([l1 1; l2 1] \ [-l1.*l1; -l2.*l2]);
    A = [0 1; -a(2) -a(1)];
    x0 = [1; 1];
    dx = @(t, x) A * x; % unstable for method 2
elseif case_num == 2
    h = 1;
    t0 = 0;
    tf = 30;
    A = [0 1; -0.025 -0.1];
    x0 = [1; 1];
    dx = @(t, x) A * x;
elseif case_num == 3
    h = 0.5;
    t0 = 0;
    tf = 30;
    A = [0 1 0; 0 0 1; -0.2 -1 -1];
    x0 = [1; 1; 1];
    dx = @(t, x) A * x;
elseif case_num == 4
    h = 0.5;
    t0 = 0;
    tf = 20;
    A = [ 998  1998; -999 -1999];
    x0 = [1; 0];
    dx = @(t, x) A * x;
    t_true = t0:1e-3:tf;
    x_true = [ 2 -1; -1  1] * exp([-1; -1000] * t_true);
elseif case_num == 5
    h = 0.5;
    t0 = 0;
    tf = 5;
    A = [ -27,  48,  -9
            5, -16,  17
          -23,  40, -17];
    x0 = [3; 0; 3];
    dx = @(t, x) A * x;
    C = [  2, 1, 2
          -1, 1, 1
           2, 1, 0];
    t_true = t0:1e-3:tf;
    x_true = C * [exp(-60.*t_true); cos(6.*t_true); sin(6.*t_true)];
elseif case_num == 6
    h = 0.01;
    t0 = 0;
    tf = 100;
    x0 = [1.2; 0];
    mu = 1e2;
    dx = @(t, x) vanderpoldemo(t, x, mu);
elseif case_num == 7
    h = 0.01;
    t0 = 0;
    tf = 11;
    x0 = [2; 0];
    mu = 1e6;
    dx = @(t, x) vanderpoldemo(t, x, mu);
elseif case_num == 8
    h = 0.001;
    t0 = 0;
    tf = 1.57;
    x0 = [1; 0];
    dx = @(t, x) -2000 .* [ cos(t .* x(1)) + sin(t .* x(2)) + 1;
                           -sin(t .* x(1)) + cos(t .* x(2)) + 1];
elseif case_num == 9
    h = 0.01;
    t0 = 0;
    % tf = 1.57;
    tf = 5;
    % x0 = [1./3; 1./2];
    % x0 = [1; 1].*1.7;
    x0 = [0.5; 0.251];
    dx = @(t, x) [(x(1).^2) .* (1./x(2) - 2) ./ (t + 1);
                  -(x(2).^2) .* (1./x(1) - 2) ./ (t + 1)];
    t_true = t0:h:tf;
    c1 = 1/x0(1) - 2;
    c2 = 2 - 1/x0(2);
    x_true = [1 ./ (c1.*cos(log(t_true+1)) + c2.*sin(log(t_true+1)) + 2);
              1 ./ (c1.*sin(log(t_true+1)) - c2.*cos(log(t_true+1)) + 2)];
elseif case_num == 10
    h = 0.01;
    t0 = 0;
    tf = 5;
    % x0 = [1./3; 1./2];
    % x0 = [1; 1].*1.7;
    x0 = [0.5; 0.2502];
    dx = @(t, x) [(x(1).^2) .* (1./x(2) - 2);
                  -(x(2).^2) .* (1./x(1) - 2)];
    t_true = t0:h:tf;
    c1 = 1/x0(1) - 2;
    c2 = 2 - 1/x0(2);
    xt_fun = @(t) [1 ./ (c1.*cos(t) + c2.*sin(t) + 2);
                   1 ./ (c1.*sin(t) - c2.*cos(t) + 2)];
    x_true = xt_fun(t_true);
elseif case_num == 11
    h = 0.01;
    t0 = 0;
    tf = 5;
    % x0 = [1./3; 1./2];
    % x0 = [1; 1].*1.7;
    % x0 = [0.501; 0.25; 0.001];
    x0 = [0.5; 0.25; 0.001];
    % x0 = [0.5; 0.2496; 0.001];
    % x0 = [1./3; 1./2; 0.01];
    a = 1e6;
    dx = @(t, x) [((x(1) - x(3)./(t + 1)).^2) .* ...
                    (1./(x(2) + x(3)./(t + 1)) - 2) ./ (t + 1) - a.*x(3);
                  -((x(2) + x(3)./(t + 1)).^2) .* ...
                    (1./(x(1) - x(3)./(t + 1)) - 2) ./ (t + 1) + a.*x(3);
                  x(3)./(t + 1) - a.*(t + 1).*x(3)];
    t_true = t0:h:tf;
    c1 = 1/(x0(1)-x0(3)) - 2;
    c2 = 2 - 1/(x0(2)+x0(3));
    c3 = a.*x0(3);
    xt_fun = @(t) [1 ./ (c1.*cos(log(t+1)) + c2.*sin(log(t+1)) + 2) + ...
                (c3./a).*exp(-a.*(0.5.*t.^2 + t));
              1 ./ (c1.*sin(log(t+1)) - c2.*cos(log(t+1)) + 2) - ...
                (c3./a).*exp(-a.*(0.5.*t.^2 + t));
              (c3./a).*(t+1).*exp(-a.*(0.5.*t.^2 + t))];
    x_true = xt_fun(t_true);
elseif case_num == 12
    h = 0.01;
    t0 = 0;
    tf = 5;
    a = 1e6;
    % x0 = [1; 1].*1.7;
    % x0 = [0.501; 0.25; 0.001];
    % x0 = [0.5; 0.2496; 0.001];
    % x0 = [1./3; 1./2; 0.01];
    x0 = [0.5; 0.25; 0.001; 0];
    % x0 = [1./3; 1./2; 0.001; 0];
    c1 = 1/(x0(1)-x0(3)) - 2;
    c2 = 2 - 1/(x0(2)+x0(3));
    c3 = a.*x0(3);
    x0(4) = -c3;
    dx = @(t, x) [((x(1) - x(3)).^2) .* ...
                    (1./(x(2) + x(3)) - 2) + x(4);
                  -((x(2) + x(3)).^2) .* ...
                    (1./(x(1) - x(3)) - 2) - x(4);
                  x(4);
                  x(4)./(sqrt(abs(1-(2./a).*log(x(3).*a./c3)))) ...
                  - a.*(sqrt(abs(1-(2./a).*log(x(3).*a./c3)))).*x(4)];
    t_true = t0:h:tf;
    xt_fun = @(t) [1 ./ (c1.*cos(t) + c2.*sin(t) + 2) + ...
                     (c3./a).*exp(-a.*(0.5.*t.^2 + t));
                   1 ./ (c1.*sin(t) - c2.*cos(t) + 2) - ...
                     (c3./a).*exp(-a.*(0.5.*t.^2 + t));
                   (c3./a).*exp(-a.*(0.5.*t.^2 + t));
                   -c3.*(t+1).*exp(-a.*(0.5.*t.^2 + t))];
    x_true = xt_fun(t_true);
elseif case_num == 13
    h = 1e4;
    t0 = 0;
    tf = 4e7;
    x0 = [1; 0; 0];
    % D = diag([1 1 1] ./ tf);
    % Jacobian is always singular!!!
    dx = @(t, y) [ (-0.04*y(1) + 1e4*y(2)*y(3))
             (0.04*y(1) - 1e4*y(2)*y(3) - 3e7*y(2)^2)
             3e7*y(2)^2 ]; % Robertson stiff problem
elseif case_num == 14
    h = 1e9;
    t0 = 0;
    tf = 1e11;
    x0 = [1; 0; 0];
    % D = diag([1 1 1] ./ tf);
    % Jacobian is always singular!!!
    dx = @(t, y) [ (-0.04*y(1) + 1e4*y(2)*y(3))
             (0.04*y(1) - 1e4*y(2)*y(3) - 3e7*y(2)^2)
             3e7*y(2)^2 ]; % Robertson stiff problem
elseif case_num == 15
    h = 1.5 / 40;
    t0 = 0;
    tf = 1.5;
    x0 = 0;
    dx = @(t,x) -2000.*(x - cos(t));
end

% % Compute exact Jacobian
% t_sym = sym('t');
% x_sym = sym('x', size(x0));
% Jf_sym = jacobian(dx(t_sym, x_sym), x_sym);
% Jf = matlabFunction(Jf_sym, 'vars', {t_sym, x_sym});

if ~exist('x_true', 'var')
    % 'true' solution
    opts = odeset('AbsTol', 1e-12, 'RelTol', 1e-8);
    [t_true, x_true] = ode15s(dx, [t0 tf], x0, opts);
    x_true = x_true.';
end

if do_comparison
    fprintf('ode23tb: ');
    tic;
    [t, x1] = ode23tb(dx, [t0 tf], x0);
    toc;
    t = t.';
    x1 = x1.';
end

if ~fixed_h; h = []; end;
fprintf('Turizo: ');
if exist('D', 'var') || fixed_h || get_hdata
    if ~exist('D', 'var')
        D = eye(length(x0));
    end
    tic;
    [t2, x2, hdata] = ercsj_scaled(dx, [t0 tf], x0, D, 1e-6, ...
        1e-3, h, get_hdata, []);
    toc;
else
    tic;
    [t2, x2] = ercsj(dx, [t0 tf], x0, 1e-6, 1e-3);
    get_hdata = false;
    toc;
end

if get_hdata
    fprintf('Total steps: %d\n', length(hdata.h_rej));
    fprintf('Rejected steps: %d\n\n', sum(hdata.h_rej));
    % figure;
    % hold all;
    % plot(hdata.ht, hdata.hcurve, 'b');
    % plot(hdata.ht(hdata.h_rej), hdata.hcurve(hdata.h_rej), 'r*');
    figure;
    hold all;
    plot(hdata.ht, 20*log(hdata.hcurve), 'b');
    plot(hdata.ht(hdata.h_rej), 20*log(hdata.hcurve(hdata.h_rej)), 'r*');
end

% Print errors
if exist('xt_fun', 'var')
    if do_comparison
        err1 = max(abs(x1 - xt_fun(t)) ./ max(abs(xt_fun(t)), 1e-3));
        err1 = max(err1(:));
        fprintf('ode23tb error: %.4f\n', err1);
    end
    err2 = max(abs(x2 - xt_fun(t2)) ./ max(abs(xt_fun(t2)), 1e-3));
    err2 = max(err2(:));
    fprintf('Turizo error: %.4f\n', err2);
end

nx = length(x0);
% nx = 0;
for i = 1:nx
    figure;
    if get_hdata
        subplot(2,1,1);
    end
    hold all;
    plot(t_true, x_true(i,:), 'b');
    if exist('x1', 'var')
        plot(t, x1(i,:), 'r-o');
        plot(t2, x2(i,:), 'g-o');
        legend({'True' 'Classic' 'New'});
    else
        plot(t2, x2(i,:), 'g-o');
        legend({'True' 'New'});
    end
    if get_hdata
        subplot(2,1,2);
        hold all;
        plot(hdata.ht,nthroot(hdata.e2(i,1:end),3), 'b');
        plot(hdata.ht,nthroot(hdata.e1(i,1:end),3), 'g');
        % plot(hdata.ht, hdata.e2(i,:), 'b');
        % plot(hdata.ht, hdata.e1(i,:), 'g');
    end
end

