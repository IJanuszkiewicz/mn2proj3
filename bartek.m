function [eigVal, eigVec, erre] = bartek(A, accuracy, maxIt, p, x0)
% Igor Januszkiewicz

if nargin < 2
    accuracy = 100*eps;
end
if nargin < 3
    maxIt = 1e4;
end
if nargin < 4
    p = 2;
end
if nargin < 5
    x0 = complex(rand(length(A), 1), rand(length(A), 1));
end

x = x0;
[L, U, P, Q] = LU(A);

x = x/norm(x, p);
i = 0;
erre = realmax;

while erre >= accuracy && i <= maxIt
    prev = x;
    x = LUsolve(L, U, P, Q, x);
    x = x/norm(x, p);
    condvec = conditionVec(x,prev);
    erre = norm(condvec, p);
    i = i + 1;
end

eigVal = x'*x/(x'*LUsolve(L, U, P, Q, x));
eigVec = x;

end % function
