function [A, lambda] = generateGoodMatrix(n, lambda1, lambda2)
% Projekt 2, zadanie 08
% Igor Januszkiewicz, 327357
%
% Funkcja tworzy macierz o zadanych najmniejszych wartościach własnych.
% W tym celu losuje reszte wartości własnych oraz macierz ortogonalną aby
% zmienić bazę

if(nargin < 2)
    lambda1 = complex(2*rand() + 1/4, 2*rand() + 1/4);
end
if(nargin < 3)
    lambda2 = lambda1*randSpin()*(1.01 + 2*rand());
end

lambdas = zeros([n,1]);
lambdas(1) = lambda1;
lambdas(2) = lambda2;
for i = 3:n
    lambdas(i) = lambda2 * randSpin() * (1.01 + 10*rand());
end

orth = randOrthMat(n);
A = orth'*diag(lambdas)*orth;
lambda = lambda1;

end % function

function c = randSpin()
c = complex(rand()-1/2, rand()-1/2);
c = c/abs(c);
end % function
