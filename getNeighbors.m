function [Y,yk] = getNeighbors(o,X,Eps1,Eps2)

d1 = sqrt((repmat(o(2),size(X,1),1)-X(:,2)).^2);
d2 = sqrt((repmat(o(1),size(X,1),1)-X(:,1)).^2);
xk = find((d1<=Eps1) & (d2<=Eps2));
yk = setdiff(xk,find(X(:,1)==o(1) & X(:,2)==o(2)));
Y = X(yk,:);
end