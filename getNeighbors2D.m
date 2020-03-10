function [Y,yk] = getNeighbors2D(o,D,Eps1,Eps2,arclen)

% d1 = sqrt((repmat(D(k,1),size(D,1),1)-D(:,1)).^2);
d1 = arclen(o(1),D(:,1))';
d2 = sqrt((repmat(o(2),size(D,1),1)-D(:,2)).^2);
xk = find((d1<=Eps1) & (d2<=Eps2));
yk = setdiff(xk,find(D(:,1)==o(1) & D(:,2)==o(2)));
Y = D(yk,:);
end