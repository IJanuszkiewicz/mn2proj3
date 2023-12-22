function [] = test4convspeed()

iters = 100;
i = 0;
maxiters = 1e4;
accuracy = 10*eps;

while i <= iters
    A = complex(rand(2), rand(2));
    [eig1, eig2] = twoDimEigFinder(A);
    if abs(eig1 - eig2) < 1000*eps()
        continue;
    end
    i = i + 1;
    eigSmaller = min(eig2, eig1, "ComparisonMethod", "abs");
    eigBigger = max(eig2, eig1, "ComparisonMethod", "abs");
    [~, ~, ~, vals] = P2Z08_IJA_OdwMetPotegowa(A, accuracy, maxiters);
    
    
    diffs = abs(eigSmaller - vals);
    logdiffs = log(diffs);
    coeff = (1:length(logdiffs))'\logdiffs;

    expectedDiv = abs(eigSmaller/eigBigger);
    actualDiv = exp(coeff);
    
    fprintf("Wartości własne: eig1 = %4f + %4fi eig2 =  %4f + %4fi\n", ...
    real(eigSmaller), imag(eigSmaller), real(eigBigger), imag(eigBigger))
    fprintf("Oczekiwany iloraz: %5f Otrzymany iloraz: %5f różnica: %e",...
        expectedDiv, actualDiv, abs(expectedDiv - actualDiv));
    fprintf("\n\n");

end  % function