function numtest1Vs(writeToFile)
% Projekt 2, zadanie 8
% Igor Januszkiewicz 327357

if nargin < 1
    writeToFile = false;
end

disp("Test porównuje różne sposoby rozwiązywania układu równań użytego w")
disp("odwrotnej metodzie potęgowej")

N = 100;

matrices = cell(N-1, 1);
startVecs = cell(N-1, 1);
lambdas = zeros(N-1, 1);
for i = 1:N - 1
    [A, lambda] = generateGoodMatrix(i + 1);
    matrices{i} = A;
    startVecs{i} = complex(rand(length(A), 1), rand(length(A), 1));
    
    lambdas(i) = lambda;
end

[time_builtin, diffs_builtin, times_builtin] = testMethod(matrices,...
    lambdas, startVecs, @builtinSolve, N);
disp("Operator A\b")
fprintf("Czas: %f, maksymalny błąd: %e\n", time_builtin,...
    max(diffs_builtin))

[time_ij, diffs_ij, times_ij] = testMethod(matrices, lambdas, startVecs,...
    @P2Z08_IJA_OdwMetPotegowa, N);
disp("Rozkład PAQ = LU")
fprintf("Czas: %f, maksymalny błąd: %e\n", time_ij, max(diffs_ij))

[time_inv, diffs_inv, times_inv] = testMethod(matrices, lambdas,...
    startVecs, @invSolve, N);
disp("Macierz odwrotna")
fprintf("Czas: %f, maksymalny błąd: %e\n", time_inv, max(diffs_inv))

[time_luandback, diffs_luandback, times_luandback] = testMethod(...
    matrices, lambdas, startVecs, @luAndBackslash, N);
disp("PAQ=LU z operatorem \")
fprintf("Czas: %f, maksymalny błąd: %e\n", time_luandback,...
    max(diffs_luandback))
if(writeToFile)
    allTimes = table((2:N)', times_builtin', times_ij', times_inv',...
        times_luandback','VariableNames', {'N', 'Operator A\b',...
        'PAQ=LU', 'Macierz odwrotna', 'PAQ=LU + operator A\b'});
    writetable(allTimes, 'solve_times.csv');
end

end % function
