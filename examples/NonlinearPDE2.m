% Nonlinear PDE:
%
% p^2 - q = y*u
%
% s.t. (x0, y0, p0, q0, u0) = (sin(s), 0, 1, 1, cos(s))

addpath('..')

f = @(x, y, p, q, u) p^2 - q - y*u;
X0 = @(s) [s, 1, s, s^2/2, s^2/2];

[x, y, p, q, u] = SolveCharacteristics(f, X0);

PlotCharacteristics(x, y, u, X0 = @(s) [s, 1, s^2/2], ...
    s_range=[-2*pi, 2*pi], t_range=[0, 5], Ns = 21, Nt = 21)

axis normal; view(-45, 45)