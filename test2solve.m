function [] = test2solve()
% Igor Januszkiewicz

n = 100;

A = rand(n);
x = rand(n,1);
b = A*x;

[L, U, P, Q] = LU(A);
resoult = LUsolve(L, U, P, Q, b);

norm(x - resoult, 2)

end % function
