function [] = test1lu()
% Projekt 2, zadanie 8
% Igor Januszkiewicz 327357

disp("Test sprawdza poprawność rozkładu PAQ = LU dla różnych macierzy")
disp("A - macierz poddana rozkładowi. Następnie zostaną wyświetlone")
disp("macierze rozkładu: najpierw policzone prze funkcje LU a potem")
disp("oczekiwany wynik.")

testCases = [
    % A, L, U, P, Q
    {[2,5;1,0], [1,0;0,1], [5,2;0,1], [1,0;0,1], [0,1;1,0]}
    {[-1,0,1;4,0,-2;1,-1,4], [1,0,0;1/4,1,0;-1/4,1/9,1], ...
    [4,-2,0;0,4.5,-1;0,0,1/9],[0,1,0;0,0,1;1,0,0], [1,0,0;0,0,1;0,1,0]}
    {[1,-2,0,1; 2,0,-1,1;-2,0,-4,2;1,1,1,1], ...
    [1,0,0,0;1/4,1,0,0;0,2/5,1,0;-1/4,2/5,-1/2,1], ...
    [-4,-2,0,2;0,2.5,0,1/2;0,0,-2,4/5;0,0,0,1.8], ...
    [0,0,1,0;0,1,0,0;1,0,0,0;0,0,0,1], ...
    [0,1,0,0;0,0,1,0;1,0,0,0;0,0,0,1]
    }
    ];

for testCase = testCases'
    A = testCase{1};
    L = testCase{2};
    U = testCase{3};
    P = testCase{4};
    Q = testCase{5};
    
    [L_calc, U_calc, P_calc, Q_calc] = LU(A, 1e-7);
    
    disp("A:")
    disp(A)
    disp("L:")
    disp(L_calc)
    disp(L)
    pause;
    disp("U:")
    disp(U_calc)
    disp(U)
    pause;
    disp("P:")
    disp(P_calc)
    disp(P)
    pause;
    disp("Q:")
    disp(Q_calc)
    disp(Q)
    disp("-------------------")
    pause;
end

end % function
