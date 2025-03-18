function X0 = GetCauchyData(f, X0, n, warn)
% Calculates p_0 and q_0 for the given set of Cauchy data for x, y, u
% X_0(s) = {x_0(s), y_0(s), u_0(s)} and outputs full Cauchy data set
% X'_0(s) = {x_0(s), y_0(s), p_0(s), q_0(s), u_0(s)}.
%
% Inputs:
% - f: A first-order PDE of the form f(x, y, p, q, u) = 0
% - X0: Initial data, paramatric curve: X0(s) = (x0, y0, u0)(s)
% - n: Number of solution for [p_0, q_0] to pick (default: 1)
% - warn: true - show warning for multiple solutions (default: true)
%
% Outputs:
% - X0: Initial data, paramatric curve: X0(s) = (x0, y0, p0, q0, u0)(s)

% Set default parameters:

if nargin < 3; n = 1; end
if nargin < 4; warn = true; end

% Define symbolic variables:

syms x y p q u s

% Convert PDE and initial data to symbolic expressions:

F(x, y, p, q, u) = f(x, y, p, q, u);

X = sym(X0(s));

x0(s) = X(1);
y0(s) = X(2);
u0(s) = X(3);

% Define system of equations for p0, q0:

syms p0 q0

Eq1 = F(x0, y0, p0, q0, u0) == 0;
Eq2 = diff(u0, s) - diff(x0, s)*p0 - diff(y0, s)*q0 == 0;

Eq = [Eq1, Eq2];

% Solve for p_0 and q_0 and warn if multiple solutions exist:

[p0, q0] = solve(Eq, p0, q0);

if length(p0) > 1 && warn

    pq_str = '';
    for i = 1:length(p0)
        pq_str = [pq_str newline 'n = ' num2str(i) ': p_0(s) = ' char(p0(i)) ', q_0(s) = ' char(q0(i))];
    end
    pq_str = [pq_str newline 'Current solution defined in X_0 is n = ' num2str(n)];

    warning(['Multiple solutions identified for p_0 and q_0. Ensure that the required solution has been identified. Solutions:' pq_str])

end

% Set p_0 and q_0 equal to specified solution:

p0(s) = p0(n);
q0(s) = q0(n);

% Output new Cauchy data X0:

X0 = matlabFunction([x0(s), y0(s), p0(s), q0(s), u0(s)]);

end