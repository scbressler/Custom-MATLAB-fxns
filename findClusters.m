function [clStat,CL] = findClusters(Xstat,thresh,df,tail)

if(strcmp(tail,'right'))
    cl = Xstat >= tinv(thresh,df);
elseif(strcmp(tail,'left'))
    cl = Xstat <= tinv(1-thresh,df);
elseif(strcmp(tail,'both'))
%     cl = abs(Xstat) >= tinv(thresh+(1-thresh)/2,df);
    cl = (Xstat<=tinv((1-thresh)/2,df) | Xstat>=tinv(thresh+(1-thresh)/2,df));
end
        
cl = [0;cl(:)]'; % zero pad for initial condition

clStat = [];
CL = [];
clNum = 0;
temp = [];

for k = 2:length(cl)
    if(cl(k))
        temp = cat(2,temp,k);
    end
    
    if(cl(k)==0 & cl(k-1)==1)
        clNum = clNum+1;
        CL{clNum} = temp-1;
        temp = [];
    end
end

if(length(CL)>0)
    for k = 1:length(CL)
        clStat(k) = sum(Xstat(CL{k}));
    end
end


% if cl(end) % Special case where last sample is a single point "cluster"
%     cl(end) = 0; % get rid of it
% end
% 
% cl1 = find(cl==1);
% cl0 = find(cl==-1)-1;
% 
% if isempty(cl1)
%         clStat = [];
%         CL = [];
%         return;
% end
% 
% %% Address special clustering cases
% if(cl0(1)<cl1(1))
%     cl0 = cl0(2:end);
% end
% 
% if(cl0(end)<cl1(end))
%     cl0(end) = length(cl);
% end
% 
% %%
% 
% 
% CL = [cl1',cl0'];
% for n = 1:size(CL,1)
%     clStat(n) = sum(Xstat(CL(n,1):CL(n,2)));
% end

end