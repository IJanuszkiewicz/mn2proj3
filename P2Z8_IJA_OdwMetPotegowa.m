function [eigVal, eigVec, erre, eigVals] = P2Z8_IJA_OdwMetPotegowa(A,...
    accuracy, maxIt, p, x0)
% Projekt 2, zadanie 8
% Igor Januszkiewicz 327357
%
% Funkcja liczy najmniejszą co do wartości bezwzględnej wartość własną
% macierzy zespolonej odwrotną metodą potęgową. Do rozwiązywania
% układów równań używany jest rozkład PAQ = LU (rozkład oparty na 
% eliminacji Gaussa z pełnym wyborem elementu głównego. 
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
%   eigVec    - Wektor własny odpowiadający obliczonej wartości własnej.
%   erre      - Oszacowanie błędu obliczania wartości własnej. Jeżeli erre
%               jest większe niż accuracy to funkcja zakończyła się w
%               skutek osiągnięcia maksymalnej liczby iteracji.
%   eigVals   - Wektor kolejnych przybliżeń wartości własnych.

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

eigVal = 1/(x'*LUsolve(L, U, P, Q, x));
eigVec = x;

end % function
