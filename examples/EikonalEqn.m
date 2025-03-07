% The Eikonal equation for a sugar pile:
%
% (p^2 + q^2 - 1) / 2 = 0
%
% (x0, y0, u0) = (a*cos(s), b*sin(s), 0)
%
% Taking (p0, q0) = (-y0'/sqrt(x0'^2 + y0'^2), x0'/sqrt(x0'^2 + y0'^2))
% gives inward pointing characteristics.

addpath('..')

f = @(x, y, p, q, u) (p^2 + q^2 - 1)/2;

a = 3; b = 2;

x0 = @(s) a*cos(s);
y0 = @(s) b*sin(s);
u0 = @(s) 0;
p0 = @(s) -b*cos(s) / sqrt(a^2 * sin(s)^2 + b^2 * cos(s)^2);
q0 = @(s) -a*sin(s) / sqrt(a^2 * sin(s)^2 + b^2 * cos(s)^2);

X0 = @(s) [x0(s), y0(s), p0(s), q0(s), u0(s)];

[x, y, p, q, u] = SolveCharacteristics(f, X0);

Ns = 21; Nt = 21; s_range = [0, 2*pi] + pi/(Ns-1);

% plot characteristics showing intersecting lines

PlotCharacteristics(x, y, u, X0 = @(s) [x0(s), y0(s), u0(s)], ...
    s_range = s_range, t_range = [0, 1.25*min(a, b)], Ns = Ns, Nt = Nt, ...
    plot_surf=false)

view(0, 90)

% Plot characteristics from the point where u = 0 to the point where xy = 0
% to ensure that the characteristics do not cross. From the first figure we
% notice that the characteristics start to cross either when x = 0 or y = 0
% so pick the first t point where that occurs via the condition xy = 0.

PlotCharacteristics(x, y, u, X0 = @(s) [x0(s), y0(s), u0(s)], ...
    s_range = s_range, Ns = Ns, Nt = Nt, ...
    t_min_cond = @(x, y, u) u, t_max_cond = @(x, y, u) x*y)

% Add top ridge where characteristics intersect

xr = @(s) max(0, (1-b^2/a^2)) * x0(s);
yr = @(s) max(0, (1-a^2/b^2)) * y0(s);
ur = @(s) min(b/a, a/b) * sqrt(a^2 * sin(s).^2 + b^2 * cos(s).^2);

hold on; fplot3(xr, yr, ur, s_range, 'k:'); view(15, 15); hold off

figure(1)

