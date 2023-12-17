function [x] = nullVec(L, U, P, Q, accuracy)
% Projekt 2, zadanie 8
% Igor Januszkiewicz 327357

x = zeros(length(U), 1);
i = length(x);

while U(i,i) < accuracy && i >= 1
    x(i) = 1;
    U(i,i) = 1;
    i = i - 1;
end

x = LUsolve(L, U, P, Q, x);

end % function