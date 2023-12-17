function [eigVal, eigVec, erre, eigVals] = P2Z08_IJA_OdwMetPotegowa(A,...
    accuracy, maxIt, p, x0)
% Projekt 2, zadanie 8
% Igor Januszkiewicz 327357
%
% Funkcja liczy najmniejszą co do wartości bezwzględnej wartość własną
% macierzy zespolonej odwrotną metodą potęgową. Do rozwiązywania
% układów równań używany jest rozkład PAQ = LU (rozkład oparty na 
% eliminacji Gaussa z pełnym wyborem elementu głównego.
%
% Parametry wejściowe:
%   A         - Kwadratowa macierz zespolona.
%   accuracy  - Dokładność wyznaczenia wartości własnej. Program zakończy
%               się kiedy wyrażenie szacujące wartość błedu wartości 
%               własnej bedzie mniejszy niż accuracy. Domyślnie wynosi
%               100*eps().
%   maxIt     - Maksymalna ilość iteracji. Domyślnie wynosi 1000
%   p         - Norma Schura, która będzie wykorzystywana do normowania
%               wektora własnego. Domyślnie wynosi 2.
%   x0        - Startowy wektor poddany iteracji. Domyślnie losowany jest
%               zespolony wektor.
% Parametry wyjściowe:
%   eigVal    - Obliczona najmniejsza co do modułu wartość własna macierzy
%               A.
%   eigVec    - Unormowany wektor własny odpowiadający obliczonej wartości
%               własnej.
%   erre      - Oszacowanie błędu obliczania wartości własnej. Jeżeli erre
%               jest większe niż accuracy to funkcja zakończyła się w
%               skutek osiągnięcia maksymalnej liczby iteracji. W przypadku
%               Macierzy osobliwej (|det(A)| < accuracy) erre wynosi
%               |det(A)|.
%   eigVals   - Wektor kolejnych przybliżeń wartości własnych.

% Ustawienie domyślnych wartośći
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

% Rozkład
[L, U, P, Q] = LU(A, accuracy);

% Sprawdzenie czy macierz A nie jest osobliwa
det = abs(prod(diag(L)*prod(diag(U))));
if det < accuracy
    eigVal = 0;
    eigVals = 0;
    erre = det;

    % Wyznaczenie wektora należącego do jądra macierzy A
    eigVec = nullVec(L, U, P, Q, accuracy);
    eigVec = eigVec/norm(eigVec, p);
    return
end

x = x/norm(x, p);
i = 2;
erre = Inf;
eigVals = zeros(maxIt, 1);
eigVals(1) = Inf;

while erre >= accuracy && i <= maxIt + 1
    prev = x;
    x = LUsolve(L, U, P, Q, x);
    eigVals(i) = 1/(prev'*x);
    x = x/norm(x, p);
    erre = abs(eigVals(i) - eigVals(i - 1));
    i = i + 1;
end

eigVals = eigVals(eigVals ~= 0);
eigVals = eigVals(2:end);
eigVal = eigVals(end);
eigVec = x;

end % function
