function [] = test3bartek()
% Igor Januszkiewicz

iters = 100;
diffs = zeros(1,iters);
accuracy = realmin;
i = 1;

while i <= iters
    A = complex(rand(2), rand(2));
    [eig1, eig2] = twoDimEigFinder(A);
    if abs(eig1 - eig2) < 1000*eps()
        continue;
    end
    eig = min(eig2, eig1);
    [eigbartek, ~, ~] = P2Z8_IJA_OdwMetPotegowa(A, accuracy);
    diffs(i) = abs(eig - eigbartek);
    i = i + 1;
end

norm(diffs, 2)

end % function
