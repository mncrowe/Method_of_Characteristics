% Damped nonlinear wave equation:
%
% u_y + u*u_x = -k*u
%
% s.t. u(x, 0) = exp(-x^2)

addpath('..')

k = 0.3;
f = @(x, y, p, q, u) q + u*p + k*u;
X0 = @(s) [s, 0, -2*s*exp(-s^2), -exp(-s^2)*(k-2*s*exp(-s^2)), exp(-s^2)];

[x, y, p, q, u] = SolveCharacteristics(f, X0);

PlotCharacteristics(x, y, u, X0 = @(s) [s, 0, exp(-s^2)], ...
    s_range=[-2, 2], t_range=[0, 2], Ns = 41, Nt = 41)

view(45, 15)