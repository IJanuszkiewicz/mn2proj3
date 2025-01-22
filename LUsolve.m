function [x] = LUsolve(L, U, P, Q, b)
% Projekt 2, zadanie 8
% Igor Januszkiewicz 327357
%
% Funkcja oblicza równanie postaci Ax = b dla przy pomocy rozkładu PAQ=LU

Pb = P * b;

n = length(b);
y = zeros(n, 1);
for i = 1:n
    y(i) = Pb(i) - L(i, 1:i-1)*y(1:i-1);
end

x_temp = zeros(n,1);
for i = n:-1:1
    x_temp(i) = (y(i) - U(i, i+1:n)*x_temp(i+1:n))/U(i, i);
end

x = Q * x_temp;

end % function
