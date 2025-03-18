% Nonlinear wave equation:
%
% u_y + u^2*u_x = 0
%
% s.t. u(x, 0) = exp(-x^2)

addpath('..')

f = @(x, y, p, q, u) q + u^2*p;
X0 = @(s) [s, 0, -2*s*exp(-s^2), 2*s*exp(-3*s^2), exp(-s^2)];

[x, y, p, q, u] = SolveCharacteristics(f, X0);

figure

syms s

Nt = 101;
t_range = linspace(0, 3, Nt);

for it = 1:Nt   % plot movie
    clf
    fplot(x(t_range(it), s), u(t_range(it), s), 'k')
    grid; xlabel('x'); ylabel('u')
    pause(0.1)
end

PlotCharacteristics(x, y, u, X0 = @(s) [s, 0, exp(-s^2)], ...
    s_range=[-2, 2], t_range=[0, 2], Ns = 41, Nt = 41)

view(15, 5); axis normal