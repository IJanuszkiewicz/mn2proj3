function [vec] = conditionVec(x,prev)

[~, index1] = max(abs(x));
[~, index2] = max(abs(prev));

vec = x*abs(x(index1))/x(index1) - prev*abs(prev(index2))/prev(index2);

end % function
