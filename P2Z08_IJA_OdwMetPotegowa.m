function [lambda, v, errEst, eigVals] = P2Z08_IJA_OdwMetPotegowa(A,...
    tol, maxIter, p, v0)
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
%             własnej bedzie mniejszy niż tol. Domyślnie wynosi
%             100*eps().
%   maxIter - Maksymalna ilość iteracji. Domyślnie wynosi 1000
%   p       - Norma Schura, która będzie wykorzystywana do normowania
%             wektora własnego. Domyślnie wynosi 2.
%   v0      - Startowy wektor poddany iteracji. Domyślnie losowany jest
%             zespolony wektor.
% Parametry wyjściowe:
%   lambda  - Obliczona najmniejsza co do modułu wartość własna macierzy
%             A.
%   v       - Unormowany wektor własny odpowiadający obliczonej wartości
%             własnej.
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
    p = 2;
end
if nargin < 5
    v0 = complex(rand(length(A), 1), rand(length(A), 1));
end

x = v0;

% Rozkład
[L, U, P, Q] = LU(A, tol);

% Sprawdzenie czy macierz A nie jest osobliwa
det = abs(prod(diag(L)*prod(diag(U))));
if det < tol
    lambda = 0;
    eigVals = 0;
    errEst = det;
    
    % Wyznaczenie wektora należącego do jądra macierzy A
    v = nullVec(L, U, P, Q, tol);
    v = v/norm(v, p);
    return
end

x = x/norm(x, p);
i = 2;
errEst = Inf;
if nargout > 3;
    eigVals = zeros(maxIter, 1);
end
eigVals(1) = Inf;

while errEst >= tol && i <= maxIter + 1
    prev = x;
    x = LUsolve(L, U, P, Q, x);
    if nargout > 3
        eigVals(i) = 1/(prev'*x);
    end
    x = x/norm(x, p);
    errEst = abs(eigVals(i) - eigVals(i - 1));
    i = i + 1;
end

eigVals = eigVals(eigVals ~= 0);
eigVals = eigVals(2:end);
lambda = eigVals(end);
v = x;

end % function
