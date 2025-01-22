function [x] = nullVec(U, Q, accuracy)
% Projekt 2, zadanie 8
% Igor Januszkiewicz 327357
%
% Funkcja oblicza wektor należący do jądra macierzy A, gdzie PAQ=LU

x = zeros(length(U), 1);
i = length(x);

while U(i,i) < accuracy && i >= 1
    x(i) = 1;
    U(i, i) = 1;
    i = i - 1;
end

x = U\x;
x = Q*x;

end % function
