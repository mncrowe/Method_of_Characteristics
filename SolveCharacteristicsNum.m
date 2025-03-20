function [x, y, p, q, u, ODE, X0] = SolveCharacteristicsNum(f, X0, s, t, nonlinear)
% Solve the given PDE, f(x, y, p, q, u) = 0, using the method of
% characteristics. Here p = u_x and q = u_y. Solutions are given
% parametrically along a set of characteristics. Characteristics start from
% the line, Gamma, prescribed by X0, hence we have imposed the condition
% that the solution surface contains this curve.
%
% Inputs:
% - f: A first-order PDE of the form f(x, y, p, q, u) = 0
% - X0: Initial (Cauchy) data along Gamma, paramatric curve of the form:
%       - X0(s) = (x0, y0, p0, q0, u0)(s)
%       - X0(s) = (x0, y0, u0)(s)
% - s: range of s to use, enter as vector of s values
% - t: range of t to use, enter as vector of t values
% - nonlinear: set to true to force method to use full Charpit equations
%       for linear problems. Is set to true for nonlinear problems.
%
% Outputs:
% - (x, y, p, q, u): values of x, y, u_x, u_y and u along characteristics
% - X0: Initial data, paramatric curve: X0(s) = (x0, y0, p0, q0, u0)(s)
%
% -------------------------------------------------------------------------
% Notes:
% The values of f and X0 should be entered as inline functions. E.g.
% F = @(x, y, p, q, u) = p - q, X0 = @(s) [s, 1, 2*s, 2*s, s^2]
%
% The initial data is required to satisfy:
%
% F(x0, y0, p0, q0, u0) = 0 and u0' = x0' * p0 + y0' * q0.
%
% This code calculates the characteristics numerically. Two methods are
% implemented here:
%
% 1. Fully nonlinear method: If the PDE is fully nonlinear or 'nonlinear'
% is set to true in the input, we use the full systme of 5 Charpit
% equations to determine the solution. This may encounter problems if there
% is a finite time blow up in p or q, i.e. wave breaking where p becomes
% infinite once the wave breaks.
%
% 2. Semilinear method: If the PDE is linear in first order derivatives of
% u, we use the 3 Monge equations for x, y, u only. The values of p and q
% are not calculated and will be output as zero arrays.
%
% -------------------------------------------------------------------------

% Define symbolic variables:

syms x y p q u r a

% Convert PDE to symbolic expression:

F(x, y, p, q, u) = f(x, y, p, q, u);

% Determine if PDE is semilinear or nonlinear:

try
    if polynomialDegree(F(x, y, p, a*p, u), p) == 1
        if nargin < 5; nonlinear = false; end
    else
        nonlinear = true;
    end
catch
    nonlinear = true;
end

% Extend Cauchy data to include p and q if only (x, y, u) is given:

if length(X0(0)) == 3
    X02 = GetCauchyData(f, X0, 1, false);
else
    X02 = X0;
end

% convert initial data to symbolic expressions:

X0s = sym(X02(r));

x0(r) = X0s(1);
y0(r) = X0s(2);
p0(r) = X0s(3);
q0(r) = X0s(4);
u0(r) = X0s(5);

% Check PDE and initial data are consistent:

if ~isAlways(F(x0, y0, p0, q0, u0) == 0)
    error('Initial data is not consistent with PDE')
end

if ~isAlways(diff(u0, r) - diff(x0, r)*p0 - diff(y0, r)*q0 == 0)
    error('Initial data in not consistent along Gamma')
end

% Define partial derivatives of F:

dFdx = diff(F, x);
dFdy = diff(F, y);
dFdp = diff(F, p);
dFdq = diff(F, q);
dFdu = diff(F, u);

% Set up Charpit's/Monge's equations and convert to a righthand-side function:

if nonlinear

    ODE = [dFdp; dFdq; -dFdx - p * dFdu; -dFdy - q * dFdu; p * dFdp + q * dFdq];
    MatFunc = matlabFunction(ODE);

    RHS = @(t, y) MatFunc(y(1), y(2), y(3), y(4), y(5));

else

    ODE(x, y, u) = [dFdp; dFdq; p * dFdp + q * dFdq - F];
    MatFunc = matlabFunction(ODE);

    RHS = @(t, y) MatFunc(y(1), y(2), y(3));

end

% Define arrays for output:

Ns = length(s);
Nt = length(t);

x = zeros(Ns, Nt);
y = zeros(Ns, Nt);
p = zeros(Ns, Nt);
q = zeros(Ns, Nt);
u = zeros(Ns, Nt);

% Solve Charpit's/Monge's equations numerically along each characteristic:

for is = 1:Ns

    [~, sol] = ode45(RHS, t, X0(s(is)));

    s_sol = size(sol); N_sol = s_sol(1);

    x(is, 1:N_sol) = sol(1:N_sol, 1);
    y(is, 1:N_sol) = sol(1:N_sol, 2);

    if nonlinear
        p(is, 1:N_sol) = sol(1:N_sol, 3);
        q(is, 1:N_sol) = sol(1:N_sol, 4);
        u(is, 1:N_sol) = sol(1:N_sol, 5);
    else
        u(is, 1:N_sol) = sol(1:N_sol, 3);
    end

end

end