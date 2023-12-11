function [x] = LUsolve(L, U, P, Q, b)
% Igor Januszkiewicz

z = L\(P*b);
y = U\z;
x = Q*y;

end % function
