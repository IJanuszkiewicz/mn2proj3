function [time, diffs, times] = testMethod(matrices, lambdas, startVecs,...
    method, N)
% Projekt 2, zadanie 8
% Igor Januszkiewicz 327357

diffs = 1:N - 1;
times = 1:N - 1;
bigTime = tic;
for i = 1:N - 1
    A = matrices{i};
    lambda = lambdas(i);
    smallTime = tic;
    [lambda_calc] = method(A, 10*eps(), 1000, startVecs{i});
    time = toc(smallTime);
    times(i) = time;
    diffs(i) = abs(lambda - lambda_calc);
end
time = toc(bigTime);

end % function
