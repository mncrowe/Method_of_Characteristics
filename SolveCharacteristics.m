function [x, y, p, q, u, F, odes, X0, n] = SolveCharacteristics(f, X0, n)
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
% - n: multiple of F to subtract from the u equation (optional)
%
% Outputs:
% - (x, y, p, q, u): values of x, y, u_x, u_y and u along characteristics
% - F: the PDE f as a symbolic expression F(x, y, p, q, u)
% - odes: Charpit's equations as symbolic expressions
% - X0: Initial data, paramatric curve: X0(s) = (x0, y0, p0, q0, u0)(s)
% - n: n value, usually highest polynomial order of terms in {p, q} in F
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
% We have subtracted n*f(x, y, p, q, u) from the final Charpit equation.
% This is valid since f = 0 and this subtraction often leads to very
% convenient cancellation for quadratic terms in p and q. The value of n is
% chosen as the polynomial order of F (in p) when q = a*p to ensure
% cancellation of highest order terms in p and q. We keep a as arbitrary to
% avoid cancellation of p terms in F. We can manually override the n value.
% -------------------------------------------------------------------------

% Extend Cauchy data to (x, y, p, q, u) if data for (x, y, u) is given:

if length(X0(0)) == 3
    X0 = GetCauchyData(f, X0, 1, false);
end

% Define symbolic variables:

syms x y p q u s

% Convert PDE and initial data to symbolic expressions:

F(x, y, p, q, u) = f(x, y, p, q, u);

X0 = sym(X0(s));

x0(s) = X0(1);
y0(s) = X0(2);
p0(s) = X0(3);
q0(s) = X0(4);
u0(s) = X0(5);

% Calculate n required to subtract n*F from the u Charpit equation:

if nargin < 3
    syms a
    try
        n = polynomialDegree(F(x, y, p, a*p, u), p);
    catch
        n = 0;
    end
end

% Check PDE and initial data are consistent:

if ~isAlways(F(x0, y0, p0, q0, u0) == 0)
    error('Initial data is not consistent with PDE')
end

if ~isAlways(diff(u0, s) - diff(x0, s)*p0 - diff(y0, s)*q0 == 0)
    error('Initial data in not consistent along Gamma')
end

% Define partial derivatives of F:

dFdx = diff(F, x);
dFdy = diff(F, y);
dFdp = diff(F, p);
dFdq = diff(F, q);
dFdu = diff(F, u);

% Set up Charpit's equations:

syms x(t) y(t) p(t) q(t) u(t)

ode1 = diff(x) == dFdp(x, y, p, q, u);
ode2 = diff(y) == dFdq(x, y, p, q, u);
ode3 = diff(p) == -dFdx(x, y, p, q, u) - p * dFdu(x, y, p, q, u);
ode4 = diff(q) == -dFdy(x, y, p, q, u) - q * dFdu(x, y, p, q, u);
ode5 = diff(u) == p * dFdp(x, y, p, q, u) + q * dFdq(x, y, p, q, u) - n*F(x, y, p, q, u);

odes = [ode1; ode2; ode3; ode4; ode5];

% Apply initial condition that characteristic starts on Gamma:

IC1 = x(0) == x0(s);
IC2 = y(0) == y0(s);
IC3 = p(0) == p0(s);
IC4 = q(0) == q0(s);
IC5 = u(0) == u0(s);

IC = [IC1, IC2, IC3, IC4, IC5];

% Solve system for characteristics, parametrised by (t, s):

syms x(t, s) y(t, s) p(t, s) q(t, s) u(t, s)

[p(t, s), q(t, s), u(t, s), x(t, s), y(t, s)] = dsolve(odes, IC);

% Simplify expressions as required:

x(t, s) = simplify(x);
y(t, s) = simplify(y);
p(t, s) = simplify(p);
q(t, s) = simplify(q);
u(t, s) = simplify(u);

end