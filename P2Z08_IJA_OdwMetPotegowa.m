function [lambda, v, errEst, eigVals] = P2Z08_IJA_OdwMetPotegowa(A,...
    tol, maxIter, v0)
% Projekt 2, zadanie 8
% Igor Januszkiewicz 327357
%
% Funkcja liczy najmniejszą co do wartości bezwzględnej wartość własną
% macierzy zespolonej odwrotną metodą potęgową. Do rozwiązywania
% układów równań używany jest rozkład PAQ = LU (rozkład oparty na
% eliminacji Gaussa z pełnym wyborem elementu głównego.
%
% Parametry wejściowe:
%   A       - Kwadratowa macierz zespolona.
%   tol     - Dokładność wyznaczenia wartości własnej. Program zakończy
%             się kiedy wyrażenie szacujące wartość błedu wartości
%             własnej bedzie mniejszy niż tol. Domyślnie wynosi 100*eps().
%   maxIter - Maksymalna ilość iteracji. Domyślnie wynosi 1000
%   v0      - Startowy wektor poddany iteracji. Domyślnie losowany jest
%             zespolony wektor.
% Parametry wyjściowe:
%   lambda  - Obliczona najmniejsza co do modułu wartość własna macierzy
%             A.
%   v       - Unormowany wektor własny odpowiadający obliczonej wartości
%             własnej
%   errEst  - Oszacowanie błędu obliczania wartości własnej. Jeżeli errEst
%             jest większe niż tol to funkcja zakończyła się w
%             skutek osiągnięcia maksymalnej liczby iteracji. W przypadku
%             Macierzy osobliwej (|det(A)| < tol) errEst wynosi
%             |det(A)|.
%   eigVals - Wektor kolejnych przybliżeń najmniejszej wartości własnej.

% Ustawienie domyślnych wartości
if nargin < 2
    tol = 100*eps;
end
if nargin < 3
    maxIter = 1e4;
end
if nargin < 4
    v0 = complex(rand(length(A), 1), rand(length(A), 1));
end

x = v0/norm(v0);

% Rozkład
[L, U, P, Q] = LU(A, tol);

% Sprawdzenie czy macierz A nie jest osobliwa
[isSingular, v, errEst] = checkNullSpace(U, Q, tol);
if(isSingular)
    lambda = 0;
    eigVals = 0;
    return
end

i = 1;
errEst = Inf;
eigVals = zeros(maxIter, 1);

while errEst >= tol && i <= maxIter
    prev = x;
    x = LUsolve(L, U, P, Q, x);
    lambda = 1/(prev'*x);
    errEst = abs(lambda)*norm(x - (1/lambda) * prev);
    x = x/norm(x);
    if nargout > 3
        eigVals(i) = lambda;
    end
    i = i + 1;
end
eigVals = eigVals(eigVals ~= 0);
v = x;
end % function
