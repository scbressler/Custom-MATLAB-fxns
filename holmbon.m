function [h,CV] = holmbon(p,a)
%
% h = holmbon(p,a)
%
% Holm-Bonferroni correction for multiple comparisons
%   INPUT VARIABLES
%       p : vector or matrix of p-values
%       a : target level of alpha (default = 0.05);
%
%   OUTPUT VARIABLE
%       h : reject null hypothesis [0:no, 1:yes]
%      CV : critical p-value cutoff; reject null hypothesis 
%           for p-values ? CV
%
% Created 2019-Nov-08: SCB
%
% Reference:
% Holm, S. (1979). A Simple Sequentially Rejective Multiple Test Procedure.
% Scandinavian Journal of Statistics, 6, 65?70.

if nargin<2
    a = 0.05; % target level of alpha (default)
end

[pI,I] = sort(p(:)); % sort p-values
N = length(p(:)); % number of p-values to consider
HB = (a./(N-(1:N)+1))'; % Holm-Bonferroni formulate
h0 = pI(:)<HB; % find all sorted p-values less than Holm-Bonferroni statistic

% Iterate through ranked p-values and compare to Holm-Bonferroni formula
% Stop at first non-rejected hypothesis

k = 0; % initialize index count
cv = 1; % intitalize while-loop state
while cv == 1
    k = k+1;
    if h0(k)
        cv = 1;
    else
        cv = 0; % first non-rejected hypothesis
    end
end

CV = pI(k); % critical p-value cutoff
h0(k:end) = 0; % set lower ranked hypothesis to non-rejected status

h(I) = h0; % re-organize hypotheses output
h = reshape(h,size(p)); % reshape to original p-value structure
