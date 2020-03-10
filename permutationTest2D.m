function [h,p,stats] = permutationTest2D(a,b,N,thresh,alpha,tail,plotIt)
%
% [h,p,stats] = permutationTest2D(a,b,N,thresh,alpha)
% Non-parametric two-dimensional (time x EEG channel location)
% cluster-based permutation testbased on:
%       Maris and Oostenveld (2007) J. Neurosci. Methods.
%
% Spatio-temporal clustering is based on ST-DBSCAN algorithm by:
%       Birant and Kut (2007) Data & Knowledge Engineering
%
% The statistical test used is the "Maximum Cluster-Level Statistic," which
% is the maximum of the cluster-level statistics, with a cluster-level
% defined as the sum of the t-statistics (single-tail, right-sided) in a
% cluster of connected samples.
%
% INPUT VARIABLES:
%       a : data set (time x channel x subject)
%       b : data set (time x channel x subject)
%           ALTERNATIVE HYPOTHESIS(H1) = a>b
%       N : number of bootstrap permutations
%  thresh : arbitrary threshold for t-statistic (default: 0.95);
%   alpha : significance level (default: 0.05)
%    tail : 'left' or 'right' tailed t-test
%
% OUTPUT VARIABLES:
%       h : hypothesis test (0=null hypothesis, 1=reject null hypothesis)
%       p : cluster-level p-value
%   stats : test summary structure
%
% Created: 2016-03-04 SCB
% Modified: 2016-05-23 SCB: added left- and right-tailed t-test capability

if nargin < 4
    thresh = 0.95; % T-distribution threshold
    alpha = 0.05;
    tail = 'right';
    plotIt = 0;
elseif nargin < 5
    alpha = 0.05;
    tail = 'right';
    plotIt = 0;
elseif nargin < 6
    tail = 'right';
    plotIt = 0;
elseif nargin < 7
    plotIt = 0;
end

%% Set test parameters
nCh = size(a,2);
nA = size(a,3);
nB = size(b,3);
df = nA-1;

if(strcmp(tail,'right'))
    THR = tinv(thresh,df);
elseif(strcmp(tail,'left'))
    THR = tinv(1-thresh,df);
end
% Cluster parameters
Eps1 = 0.5;     % Maximum spatial distance for cluster inclusion
Eps2 = 1;       % Maximum temporal distance for cluster inclusion
MinPts = 1;     % Minimum number of points within Eps1 and Eps2 distance

%% Collect channel location data (default: 32 channels)
switch nCh
    case 32
        load /Users/scbressler/Documents/MATLAB/chanlocs32.mat;
    case 64
        load /Users/scbressler/Documents/MATLAB/chanlocs64.mat;
end

lat = [chanlocs.sph_theta];
lon = [chanlocs.sph_phi];
try
    SPH = referenceSphere('unit sphere');

    % Calculate arc length distances between channels
    for m = 1:nCh
        for n = 1:nCh
            [arclen(m,n),az(m,n)] = distance(lat(m),lon(m),lat(n),lon(n),SPH);
        end
    end
catch
    switch nCh
    case 32
        load /Users/scbressler/Documents/MATLAB/customfxns/UnitSphere_32ch.mat;
    case 64
        load /Users/scbressler/Documents/MATLAB/customfxns/UnitSphere_64ch.mat;
    end
end

% Calculate the OBSERVED test statistic and relative channel distances
for n = 1:nCh
    [~,~,~,statsOb] = ttest(squeeze(a(:,n,:))',squeeze(b(:,n,:))','Tail',tail);
    Xstat(n,:) = statsOb.tstat;
end

[r,c] = find(double(Xstat>=THR));
D = [r,c];
nD = size(D,1);
CL = zeros(nD,1);

%% Get observed clusters
[ObsClstrStat,ObsClstr] = findClusters2D(Xstat,arclen,thresh,df,tail,Eps1,Eps2,MinPts);

if(isempty(ObsClstrStat));
    h = zeros(1,size(a,2));
    Ci = [];
    p = [];
    stats = struct('ClusterStats',[],'ClusterIndices',Ci,'pValue',[],...
                   'p',p,'threshold',thresh,'alpha',alpha,...
                   'numPermutations',N,'method',sprintf('t-Tail:%s',tail));
    fprintf('ERROR: No clusters in OBSERVED DATA. Adjust threshold variable\n');
    return;
end

%% Permutation step
poolobj = parpool;

% Concatenate the data
Xcat = cat(3,a,b);

parfor k = 1:N
    % Permute the data by subjects
    Ik = randperm(nA+nB);
    Ystat = [];
    
    for n = 1:nCh
        [~,~,~,Xtest] = ttest(squeeze(Xcat(:,n,Ik(1:nA)))',...
            squeeze(Xcat(:,n,Ik(nA+1:end)))',...
            'Tail',tail);
%         Ystat(n,:) = Xtest.tstat;
        Ystat = cat(1,Ystat,Xtest.tstat);
        
    end
        
    
    [clstr{k},kCL] = findClusters2D(Ystat,arclen,thresh,df,tail);
      % Delete later if below works
%     try prmStat(k) = max(clstr{k});
%     catch prmStat(k) = NaN;
%     end
    
    try 
        if(strcmp(tail,'left'))
            prmStat(k) = min(clstr{k});
        elseif(strcmp(tail,'right'));
            prmStat(k) = max(clstr{k});
        end
    catch prmStat(k) = NaN;
    end
end
delete(poolobj);

%% Test for significance based on null distribution derived from permuted data
h = zeros(size(Xstat));
[f,x] = ecdf(prmStat);
pVal = x(find(f>=(1-alpha),1,'first'));

for n = find(ObsClstrStat>=pVal)
    h(sub2ind(size(h),ObsClstr{n}(:,1),ObsClstr{n}(:,2))) = 1;
end

for n = 1:length(ObsClstrStat)
    if(ObsClstrStat(n)>x(end))
        p(n) = 0;
    else
        p(n) = 1-f(find(x>=ObsClstrStat(n),1,'first'));
    end
end

ClstrMatrix = zeros(size(Xstat));
for k = 1:length(ObsClstr)
    for j=1:size(ObsClstr{k},1)
        point=ObsClstr{k}(j,:);
        ClstrMatrix(point(1),point(2))=k;
    end
end

stats.ClusterStats = ObsClstrStat;
stats.ClusterIndices = ObsClstr;
stats.ClusterMatrix = ClstrMatrix;
stats.pValue = pVal;
stats.p = p;
stats.threshold = thresh;
stats.alpha = alpha;
stats.numPermuations = N;
stats.method = 'One-Tailed:right';


%% Plots
if plotIt
    [T,CH] = meshgrid(1:200,1:nCh);
    
    figure(1);
    s1 = surf(T,CH,Xstat);
    hold on;
    s2 = surf(T,CH,repmat(tinv(0.95,12),nCh,200),'LineStyle','none','FaceAlpha',0.8);
    s2.CData = -6*ones(nCh,200);
    xlabel('Sample'); ylabel('Channel#'); zlabel('t-score');
    
    figure(2);
    set(2,'Position',[-1599 1 1600 928]);
    
    colormap(colorcube);
    subplot(1,2,1);
    pcolor(double(Xstat>=THR));
    xlabel('Time(sec)'); ylabel('Channel #');
    
    subplot(1,2,2);
    pcolor(ClstrMatrix); colorbar;
    xlabel('Sample'); ylabel('Channel #');
    
    figure(3);
    set(3,'Position',[83 365 1081 420]);
    pause;
    for K = 1:size(Xstat,2)
        clf
        subplot(1,2,1);
        topoplot(mean((a(K,:,:)-b(K,:,:)),3).*h(:,K)',...
                 chanlocs,'electrodes','numbers',...
                 'maplimits',0.5*[-1 1]);
        title(sprintf('Sample:%1.0f',K));
        colorbar;
        
        subplot(1,2,2);
        hNaN = h;
        hNaN(h==0) = NaN;
        pc3 = pcolor(hNaN.*ClstrMatrix);
        set(pc3,'LineStyle','none');
        xlabel('Time(sec)'); ylabel('Channel #');
        line(K*[1 1],ylim,'Color','c','LineWidth',2,'LineStyle','--');
        pause(0.05);
    end
    
    pause(1);
end
        
end


    