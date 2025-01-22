function [] = test4convspeed()
% Projekt 2, zadanie 8
% Igor Januszkiewicz 327357

disp("Test sprawdza szybkość zbieżności funkcji. Błąd obliczenia wartości")
disp("powinien zbiegać do 0 jak |lambda_1/lambda_2|^k, gdzie lambda_1 to")
disp("najmniejsza co do wartości bezwzględnej wartość własna, lambda_2")
disp("to druga wartość własna. iloraz znajdywany jest przy pomocy")
disp("dopasowania liniowego logarytmu z kolejnych błędów obliczeń.")
fprintf('\n')
iters = 4;
i = 0;
maxiters = 1e3;
accuracy = 1e-9;

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
    
    fprintf("Wartości własne: eig1=%.4f+%.4fi eig2=%.4f+%.4fi\n", ...
        real(eigSmaller), imag(eigSmaller), real(eigBigger),...
        imag(eigBigger))
    fprintf("Oczekiwany iloraz: %5f Otrzymany: %5f różnica: %e",...
        expectedDiv, actualDiv, abs(expectedDiv - actualDiv));
    fprintf("\n\n");
    
end  % function
