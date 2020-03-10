function [h,p,stats] = permutationTest(a,b,N,thresh,alpha,tail)
%
% [h,p,stats] = permutationTest(a,b,N,thresh,alpha,tail)
%   Performs cluster-based nonparametric statistics on two times series
%   across two different conditions.  Assumes same subjects in both
%   conditions.  The null distribution of the test statistic is derived
%   from N bootstrapped permutations of the data.
%
%   The test statistic used in this version is the sum of the t-statistic
%   within a time cluster that exceeds a certain thresholds, thresh,
%   defined by the user.  The choice of the threshold value is arbitrary,
%   but does have implications in terms of the sensitivity of the test.
%   It does not affect False Alarm rate of the test.
%
%   For a more detailed description of the methods and statistical validity
%   see: 
%
%   Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing 
%       of EEG- and MEG-data. Journal of Neuroscience Methods, 164(1), 
%       177?190. http://doi.org/10.1016/j.jneumeth.2007.03.024
%
% INPUT VARIABLES
%       a : multi-subject time series for Condition A <subj x samples>
%       b : multi-subject time series for Condition B <subj x samples>
%       N : number of bootstrap permutations (recommend at least 1000)
%  thresh : user selected t-statistic threshold (default=0.95)
%   alpha : p-value criteria
%    tail : direction of the one-sided t-test ['left' or 'right']
%
% OUTPUT VARIABLE
%       h : sample-by-sample test decision variable 
%               1=reject null hypothesis at significance level, alpha
%               0=failure to reject null hypothesis
%       p : cluster-by-cluster p-values
%   stats : structure of test results
%
% Created by: Scott Bressler 08-Jun-2016
% Updated: 21-May-2019, now can do left, right, and two-tailed t-test and
%                       can now do paired and two-sample t-tests

if nargin < 4
    thresh = 0.95; % t-distribution threshold
    alpha = 0.05;  % significance level
    tail = 'left'; % t-test tail
elseif nargin < 5
    alpha = 0.05;
    tail = 'left';
elseif nargin < 6
    tail = 'left';
end

% Initialize test decision variable, h
h = zeros(1,size(a,2));

% Determine number of subjects per condition
nA = size(a,1);
nB = size(b,1);
df = nA+nB-2;  % degrees of freedom

% Determine correct t-stat threshold based on tail
if(strcmp(tail,'right'))
    THR = tinv(thresh,df);
elseif(strcmp(tail,'left'))
    THR = tinv(1-thresh,df);
end

% Id = cat(1,ones(nA,1),2*ones(nB,1)); % delete later

% Calculate the "Observed" test statistic
if(nA==nB) % use paired t-test
    fprintf('Running paired t-test...');
    [Xobs.h,Xobs.p,Xobs.ci,Xobs.stats] = ttest(a,b,'Tail',tail);
    df = Xobs.stats.df(1);
else % use two-sample t-test
    fprintf('Running two-sample t-test...');
    [Xobs.h,Xobs.p,Xobs.ci,Xobs.stats] = ttest2(a,b,'Tail',tail);
    df = Xobs.stats.df(1);
end
[ObsClstrStat,kClstr] = findClusters(Xobs.stats.tstat,thresh,df,tail);
fprintf('DONE\n');

% Return null results if no clusters in Observed Datac
if(isempty(ObsClstrStat));
    h = zeros(1,size(a,2));
    Ci = [];
    p = [];
    stats = struct('ClusterStats',[],'ClusterIndices',Ci,'pValue',[],...
                   'p',p,'threshold',thresh,'alpha',alpha,...
                   'numPermutations',N,'method',sprintf('t-Tailed:%s',tail));
    fprintf('WARNING: No clusters in OBSERVED DATA.\n');
    fprintf('         Consider adjusting threshold variable\n');
    return;
end

%% Prepare data for bootstrap permutation procedure
Xcat = cat(1,a,b); % gather both groups into a single matrix

Xtest = repmat([],N,size(a,2)); % initialize variable

for k = 1:N % loop through all bootstrap trials
    Ik = randperm(nA+nB); % randomly permute subject data
    
    % Calculate permutation test statistic (difference a-b)
    if(nA==nB) % paired t-test
        [~,~,~,Xtest] = ttest(Xcat(Ik(1:nA),:),Xcat(Ik(nA+1:end),:),'Tail',tail);
    else % two-sample t-test
        [~,~,~,Xtest] = ttest2(Xcat(Ik(1:nA),:),Xcat(Ik(nA+1:end),:),'Tail',tail);
    end
    % Locate time serires clusters
    [clstr,kCL] = findClusters(Xtest.tstat,thresh,df,tail);
    
    % Determine max/min cluster-level permuted test statistic
    try 
        if(strcmp(tail,'left'))
            prmStat(k) = min(clstr);
        elseif(strcmp(tail,'right'))
            prmStat(k) = max(clstr);
        elseif(strcmp(tail,'both'))
            [~,I] = max(abs(clstr));
            prmStat(k) = clstr(I);
        end
    catch prmStat(k) = NaN;
    end
    
end

% Calculate permuation null distribution, based on empirical cumulative
% distribution function
if(strcmp(tail,'both'))
    [f,x] = ecdf(abs(prmStat));
else
    [f,x] = ecdf(prmStat);
end

% Determine alpha significance value, pVal, and test decision variable
switch tail
    case 'right'
        tVal = x(find(f>=(1-alpha),1,'first'));
        if(sum(ObsClstrStat>=tVal))
            for n = find(ObsClstrStat>=tVal)
                h(kClstr{n}) = 1;
            end
        end
    case 'left'
        tVal = x(find(f<=alpha,1,'last'));
        if(sum(ObsClstrStat<=tVal))
            for n = find(ObsClstrStat<=tVal)
                h(kClstr{n}) = 1;
            end
        end
    case 'both'
        tVal = x(find(f>=(1-alpha)+(alpha/2),1,'first'));
%         tVal(2) = x(find(f<=(alpha/2),1,'last'));
        if(sum(abs(ObsClstrStat)>=tVal))
            for n = find(abs(ObsClstrStat)>=tVal)
                h(kClstr{n}) = 1;
            end
        end
end


% Determine cluster-level p-values
for n = 1:length(ObsClstrStat)
    switch tail
        case 'right'
            if(ObsClstrStat(n)<min(x))
                p(n) = 0;
            elseif(ObsClstrStat(n)>max(x))
                p(n) = 1;
            else
                p(n) = f(find(x>=ObsClstrStat(n),1,'first'));
            end
        case 'left'
            if(ObsClstrStat(n)>max(x))
                p(n) = 1;
            elseif(ObsClstrStat(n)<min(x))
                p(n) = 0;
            else
                p(n) = 1-f(find(x>=ObsClstrStat(n),1,'first'));
            end
        case 'both'
            if(abs(ObsClstrStat(n))<min(x))
                p(n) = 1;
            elseif(abs(ObsClstrStat(n))>max(x))
                p(n) = 0;
            else
                p(n) = 1-f(find(x>=abs(ObsClstrStat(n)),1,'first'));
            end
    end
end

% Output stats summary structure
stats.ClusterStats = ObsClstrStat;
stats.ClusterIndices = kClstr;
stats.tValue = tVal;
stats.p = p;
stats.df = df;
stats.threshold = thresh;
stats.alpha = alpha;
stats.numPermuations = N;
stats.method = sprintf('t-Tailed:%s',tail);
        
end


    