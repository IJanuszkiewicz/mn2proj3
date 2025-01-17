function [L, U, P, Q] = LU(A, tol)
% Projekt 2, zadanie 08
% Igor Januszkiewicz, 327357
%
% Funkcja oblicza rozkład PAQ = LU.

U = A;
P = eye(size(A));
Q = P;
L = zeros(size(U));

for k = 1:length(U) - 1
    % wybór
    [maxEl, i] = max(abs(U(k:end, k:end)), [], "all");
    if maxEl < tol
        L = L + eye(size(L));
        return
    end
    [i, j] = ind2sub(size(U(k:end, k:end)), i);
    i = i + k - 1;
    j = j + k - 1;
    
    % zamiana
    U(:, [k, j]) = U(:, [j, k]);
    Q(:, [k, j]) = Q(:, [j, k]);
    U([k, i], :) = U([i, k], :);
    P([k, i], :) = P([i, k], :);
    L([k, i], :) = L([i, k], :);
    
    % odejmowanie
    for l = k + 1:length(U)
        coeff = U(l, k)/U(k, k);
        L(l, k) = coeff;
        U(l, k:end) = U(l, k:end) - coeff*U(k, k:end);
    end
    U(k + 1:end, k) = 0;
end

L = L + eye(size(L));

end % function
