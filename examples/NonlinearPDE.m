% Nonlinear PDE:
%
% p*q + log(x) = 0
%
% s.t. (x0, y0, p0, q0, u0) = (x0, s, -log(x0), 1, s)

addpath('..')

x0 = 1;
f = @(x, y, p, q, u) p*q + log(x);
X0 = @(s) [x0, s, -log(x0), 1, s];

[x, y, p, q, u] = SolveCharacteristics(f, X0);

PlotCharacteristics(x, y, u, X0 = @(s) [x0, s, s], ...
    s_range=[-2, 2], t_range=[0, 2], Ns = 21, Nt = 21)

axis normal; view(-45, 45)
