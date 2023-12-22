function [A] = randOrthMat(n, tol)

if nargin < 2
    tol = 500*eps;
end
A = zeros(n);

randvec = @() complex(rand(n,1)*2 - 1, rand(n,1)*2 - 1);

vi = randvec();
A(:,1) = vi./norm(vi);

for i = 2:n
    vinorm = 0;
    while vinorm < tol
        vi = randvec();
        vi = vi - A(:, 1:i-1)*(A(:,1:i-1)' *vi);
        vinorm = norm(vi);
    end
    A(:, i) = vi./vinorm;
end

end % function