function [eigval1, eigval2] = twoDimEigFinder(A)
% Igor Januszkiewicz

a = A(1,1);
b = A(1,2);
c = A(2,1);
d = A(2,2);
eigval1 = (a + d - sqrt((a + d)*(a + d) - 4*(a*d - c*b)))/2;
eigval2 = (a + d + sqrt((a + d)*(a + d) - 4*(a*d - c*b)))/2;

end % function
