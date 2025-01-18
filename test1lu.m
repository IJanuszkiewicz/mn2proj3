function [] = test1lu()
% Igor Januszkiewicz
text = "Test sprawdza poprawność rozkładu PAQ = LU dla różnych macierzy"

testCases = [
    % A, L, U, P, Q
    {[2,5; 1,0], [1,0; 0,1], [2,5; 0,1], [1,0; 0,1], [0,1; 1,0]}
    ];

for testCase = testCases'
    A = testCase{1};
    L = testCase{2};
    U = testCase{3};
    P = testCase{4};
    Q = testCase{5};
    
    [L_calc, U_cals, P_calc, Q_calc] = LU(A, 1e-7);
    
    disp(A)
    
end

end % function
