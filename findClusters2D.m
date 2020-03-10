function [clStat,C] = findClusters2D(Xstat,arclen,thresh,df,tail,Eps1,Eps2,MinPts)

if nargin<6
    Eps1 = 0.5;   % Maximum spatial distance value
    Eps2 = 1;   % Maximum temporal distance value
    MinPts = 1; % Minimum number of point within Eps1 and Eps2 distance
end

%%
if(strcmp(tail,'right'))
    [r,c] = find(double(Xstat>=tinv(thresh,df)));
elseif(strcmp(tail,'left'))
    [r,c] = find(double(Xstat<=tinv(1-thresh,df)));
end

D = [r,c];
nD = size(D,1);
CL = zeros(nD,1);

% Initialize parameters
clN = 0; % Cluster number

for k = 1:nD
    if(CL(k)==0)
        % Retrieve nearest neighbors, X
        [X,xk] = getNeighbors2D(D(k,:),D,Eps1,Eps2,arclen);
          
        if size(X,1)<MinPts % IF number of neigbors<MinPts
            CL(k) = 0; % Mark as noise
        else
            clN = clN+1; % Define new cluster
           % Mark nearest neighbors with current CLUSTER label "clN"
            CL(xk) = clN;
           % Create a 'stack' object
            stack=X;
          
           while ~isempty(stack)
                    currentObj=stack(1,:); 
                    %need to take this out of the stack, equivalent to
                    %pop()
                    stack(1,:)=[];
                    % Iterate through stack and find nearest neighbors
                    [Y,yk] = getNeighbors2D(currentObj,D,Eps1,Eps2,arclen);
                    
                    if(size(Y,1)>=MinPts)
                       for p=1:size(Y,1)
                           if CL(yk(p))==0
                              CL(yk(p))=clN;
%                        CL(yk(CL(yk)==0)) = clN;
                              stack=[Y(p,:); stack]; %push
                           end
                       end
                    end
           end
        end
    end
end

if(unique(CL)==0) % if no clusters found, send empty variable
    C = [];
    clStat = [];
else
    for k = 1:length(unique(CL))
        C{k} = D(CL==k,:);
        clStat(k) = sum(Xstat(sub2ind(size(Xstat),C{k}(:,1),C{k}(:,2))));
    end
end


