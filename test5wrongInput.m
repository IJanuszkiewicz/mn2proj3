function [] = test5wrongInput()
% Projekt 2, zadanie 8
% Igor Januszkiewicz 327357

disp("Test sprawdza zachowanie funkcji w przypadku macierzy nie")
disp("spełniających warunków zbieżności.")
disp(' ')

disp("Macierz osobliwa:")
A = [1,2,3; -1,2,1; 2,4,6];
disp(A)
[lambda] = P2Z08_IJA_OdwMetPotegowa(A);
fprintf("Obliczona wartość własna: %f\n", lambda);

disp("Macierz bez najmniejszej wartości własnej:");
A = [2,1,0; 0,-2,1; 0,0,sqrt(2) + 1i*sqrt(2)];
disp(A)
[lambda] = P2Z08_IJA_OdwMetPotegowa(A);
fprintf("Obliczona wartość własna: %f\n", lambda);

end % function
