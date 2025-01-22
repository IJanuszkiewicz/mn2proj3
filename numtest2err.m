function [] = numtest2err(writeToFile)
% Projekt 2, zadanie 8
% Igor Januszkiewicz 327357

disp("Funkcja sprawdza zależność czasu obliczeń od wielkości macierzy i ")
disp("zadanej dokładności")

if(nargin < 1)
    writeToFile = false;
end

Ns = [10, 100, 200];

numTols = 20;
minTol = -18;
maxTol = -12;
tols = logspace(minTol, maxTol , numTols);

times = zeros(length(Ns), length(tols));

for n = 1:length(Ns)
    matrix = generateGoodMatrix(Ns(n));
    for tol = 1:length(tols)
        tic;
        P2Z08_IJA_OdwMetPotegowa(matrix, tols(tol));
        times(n, tol) = toc;
    end
end

loglog(tols, times);
title("Czas obliczeń a zadany błąd")
xlabel("tolerancja")
ylabel("czas obliczeń [s]")


if(writeToFile)
    resultTable = table(tols', times(1,:)', times(2,:)', times(3,:)',...
        'VariableNames', {'tolerances', 'n1', 'n2', 'n3'});
    writetable(resultTable, 'errors.csv');
end

end % function
