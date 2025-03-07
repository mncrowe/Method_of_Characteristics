function PlotCharacteristics(x, y, u, options)
% Plots characteristics (x, y, u) in 3D space
%
% Inputs:
% - (x, y, u): characteristics as symbolic expressions of t and s
% - options;
%   - s_range: range of parameter s to plot for (default: [-1, 1])
%   - t_range: range of parameter t to plot for (default: [0, 1])
%   - Ns: number of characteristics to plot (default: 21)
%   - Nt: gridpoint of t for surface plot (default: 101)
%   - X0: (x0, y0, u0)(s) as an inline function of s, optional (default: 0)
%   - plot_chars: if 'true', plots characteristics (default: true)
%   - plot_surf: if 'true', plots surface by interpolation (default: true)
%   - (t_min_cond, t_max_cond): conditions on (x,y,p,q,u) for t_range
%
% Note: this function is designed to plot solutions from the output of
% 'SolveCharacteristics.m' so the symbolic inputs are consistent with the
% outputs from that function.

% Define input arguments:

arguments
    x
    y
    u
    options.s_range = [-1 1]
    options.t_range = [0 1]
    options.Ns = 21
    options.Nt = 101
    options.X0 = 0
    options.plot_chars = true
    options.plot_surf = true
    options.t_min_cond = 0
    options.t_max_cond = 0
end

% Define internal parameters:

syms t s

s_range = linspace(options.s_range(1), options.s_range(2), options.Ns);

% Create figure:

figure
hold on

% Define Gamma if specified as an option:

if ~isequal(options.X0, 0)

    X0 = sym(options.X0(s));

    x0(s) = X0(1);
    y0(s) = X0(2);
    u0(s) = X0(3);

    fplot3(x0(s), y0(s), u0(s), [s_range(1), s_range(end)], 'k--')

end

% Parse conditions for start and end of characteristics:

if ~isequal(options.t_min_cond, 0)
    syms xs ys us
    t_min_cond(xs, ys, us) = options.t_min_cond(xs, ys, us);
end

if ~isequal(options.t_max_cond, 0)
    syms xs ys us
    t_max_cond(xs, ys, us) = options.t_max_cond(xs, ys, us);
end

% Express x, y, u as functions of t, s

x(t, s) = x;
y(t, s) = y;
u(t, s) = u;

% Calculate start and end of each characteristic:

t_range = zeros(options.Ns, 2);

for is = 1:options.Ns
    
    if ~isequal(options.t_min_cond, 0)
        t_min = max(solve(t_min_cond(x(t, s_range(is)), y(t, s_range(is)), u(t, s_range(is)))));
    else
        t_min = options.t_range(1);
    end
    
    if ~isequal(options.t_max_cond, 0)
        t_max = min(solve(t_max_cond(x(t, s_range(is)), y(t, s_range(is)), u(t, s_range(is)))));
    else
        t_max = options.t_range(2);
    end

        t_range(is, :) = double([t_min, t_max]);

end

% Plot surface by interpolation over characteristics

if options.plot_surf

    xM = zeros(options.Ns, options.Nt);
    yM = zeros(options.Ns, options.Nt);
    uM = zeros(options.Ns, options.Nt);
    
    for is = 1:options.Ns

        t_points = linspace(t_range(is, 1), t_range(is, 2), options.Nt);

        xM(is, :) = x(t_points, s_range(is));
        yM(is, :) = y(t_points, s_range(is));
        uM(is, :) = u(t_points, s_range(is));

    end
    
    surf(xM, yM, uM); shading interp

end

% Plot characteristics

if options.plot_chars

    for is = 1:options.Ns

        fplot3(x(t, s_range(is)), y(t, s_range(is)), u(t, s_range(is)), t_range(is, :), 'k');

    end

end

% Set axis labels and rescale figure

xlabel('x'); ylabel('y'); zlabel('u'); axis equal; grid

hold off

end