function [isSingular, v, errEst] = checkNullSpace(U,Q,tol)
% Projekt 2, zadanie 8
% Igor Januszkiewicz 327357
%
% funkcja sprawdza czy macierz nie jest osobliwa, jeśli jest to zwraca
% wektor z jądra tej macierzy

isSingular = false;
v = [];
errEst = inf;
det = abs(prod(diag(U)));
if det < tol
    errEst = det;
    isSingular = true;
    
    % Wyznaczenie wektora należącego do jądra macierzy A
    v = nullVec(U, Q, tol);
    v = v/norm(v);
    return
end

end % function
