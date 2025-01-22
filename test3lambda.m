function [] = test3lambda()
% Projekt 2, zadanie 8
% Igor Januszkiewicz 327357

tol = 100*eps();
disp("Test sprawdza poprawność obliczania najmniejszej co do modułu")
disp("wartości własnej. Zostanie wyświetlona analizowana macierz A, jej")
disp("najmniejsza co do modułu wartość własna lambda, błąd obliczeń")
disp("err = |lambda - lambda_obliczona|");
testCases = [
    % A, lambda
    {[2,1;1,2], 1}
    {[2,0,0; 0,3,4;0,4,9], 1}
    {[1,0,0;1,2,0;2,3,3], 1}
    {[3,1,2;2,4,4;1,1,6], 2}
    ];


i = 0;
for testCase = testCases'
    A = testCase{1};
    lambda = testCase{2};
    [lambda_calc] = P2Z08_IJA_OdwMetPotegowa(A, tol);
    i = i + 1;
    if i == 3
        pause;
    end
    
    disp("A: ")
    disp(A)
    err = abs(lambda - lambda_calc);
    fprintf("lambda = %f, lambda_obliczona = %f, err = %e\n", lambda,...
        lambda_calc, err);
    disp("---------------------------------------------------------------")
end

end % function
