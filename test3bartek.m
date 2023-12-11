function [] = test3bartek()
% Igor Januszkiewicz

iters = 10;
diffs = zeros(1,iters);
accuracy = realmin;
stopByCond = diffs;
i = 1;

while i <= iters
    A = complex(rand(2), rand(2));
    [eig1, eig2] = twoDimEigFinder(A);
    if abs(eig1 - eig2) < 1000*eps()
        continue;
    end
    eig = min(eig2, eig1);
    [eigbartek, ~, ~] = bartek(A, accuracy);
    diffs(i) = abs(eig - eigbartek);
    i = i + 1;
end

norm(diffs, 2)

end % function
