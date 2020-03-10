function h = jitterPlot(X,Y,Color,MarkerType,jitX)
% h = jitterPlot(X,MarkerType,Color)
%   Used to overlay on top of box plots to show individual subject data.
%   Duplicate data points are jittered horizontally based on a fixed unit
%   (dx)
%
% INPUT VARIABLES
%            X : position on the x-axis
%            Y : data organized by columns
%        Color : 3-element vector of color
%   MarkerType : string defining marker type (e.g. 'o')
%         jitX : amount of jitter for duplicates
%
% OUTPUT VARIABLE
%            h : plot handle
%
% Created 2018-Sep-05 SCB

if nargin < 3
    Color = [0.6 0.6 0.6];
    MarkerType = 'o';
    jitX = 0.01;
elseif nargin < 4
    MarkerType = 'o';
    jitX = 0.01;
elseif nargin < 5
    jitX = 0.01;
end

for k = 1:size(Y,2)
    lvl = unique(Y(:,k));
    lvl = lvl(~isnan(lvl));
    
    for m = 1:length(lvl)
        N = sum(Y(:,k)==lvl(m)); %number of duplicate values
        dx = jitX*((0:N-1)-(N-1)/2); % X-axis offset
        h(k) = plot(dx+X(k),Y(Y(:,k)==lvl(m),k),MarkerType); hold on;
        h(k).MarkerSize = 4;
        h(k).Color = Color;
    end
end
