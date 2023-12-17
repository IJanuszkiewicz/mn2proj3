function [] = numtestErr()

n = 100;
iters = 10;

for i = 1:iters
    A = complex(rand(n)*10, rand(n)*10);
    [~, ~, err] = bartek(A);
    err
end

end % function


% zrobić testy dla macierzy clementa wilkinsona!!