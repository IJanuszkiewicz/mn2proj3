function [] = test2solve()
% Projekt 2, zadanie 8
% Igor Januszkiewicz 327357

disp("Test sprawdzą poprawność rozwiązywania równań Ax = b przy pomocy")
disp("rozkładu PAQ = LU. Zostanie wyświetlona macierz A, wektor x,")
disp("rozwiązanie równania przy pomocy funkcji LUsolve x_calculated oraz")
disp("wartość błędu err = ||x - x_calculated||")

testCases = [
    % A, x
    {[1,2;-2,5], [1;-1]}
    {[1,-1,5;2,3,4;-3,4,1], [1;2;3]}
    {[complex(5,-2),complex(0,2);complex(4,-5),complex(-1,2)],...
    [complex(0,2);complex(-2,1)]}
    ];

for testCase = testCases'
    A = testCase{1};
    x = testCase{2};
    
    b = A*x;
    [L, U, P, Q] = LU(A, 1e-10);
    result = LUsolve(L, U, P, Q, b);
    disp("A:")
    disp(A)
    disp("x:")
    disp(x)
    disp("x_calculated:")
    disp(result)
    disp("err:")
    disp(norm(x - result, 2))
    disp("--------------------")
    pause;
end

end % function
