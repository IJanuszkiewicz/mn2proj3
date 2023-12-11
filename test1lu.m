function [] = test1lu()
% Igor Januszkiewicz

A = magic(4);
As = sparse(A);
  [L, U, P, Q] = LU(A);
[Ls, Us, Ps, Qs] = lu(As);

Ls = full(Ls)
L
Us = full(Us)
U
Ps = full(Ps)
P
Qs = full(Qs)
Q

P*A*Q - L*U


end % function
