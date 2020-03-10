function w = aespa(x,y)
%
% w = aespa(x,y)
%
% Auditory-Evoked Spread Spectrum Analysis algorithm (AESPA): estimates the
% "impulse response" of the auditory system from a known acoustic
% stimulus, x, and the EEG-measured scalp voltages, y, using linear
% least-squares regression.
%
% Based on the model:
%   y(t) = w(tau) * x(t) + noise
%
% Solution for w(tau) is based on the analytic solution:
%   w = <xx'><xy>, where
%       <xx'> is the autocorrelation matrix of x
%       <xy> is the cross correlation of x and y
%
% INPUT VARIABLES
%   x : regressor, acoustic input signal [N samples x M trials] matrix
%   y : response, EEG channel scalp voltage [N samples x M trials] matrix
%
% OUTPUT VARIABLE
%   w : estimated impulse response function [N samples x M trials] matrix
%
% References:
% RK Maddox and AKC Lee (2018) "Auditory Brainstem Responses to Continuous
% Natural Speech in Human Listeners" eNeuro, 5(1), e0441-17.2018 1?13
%
% EC Lalor, et al., (2009) "Resolving Precise Temporal Processing
% Properties of the Auditory System Using Continuous Stimuli" J
% Neurophysiol, 102: 349-359.
%
% EC Lalor, et al., (2006) "The VESPA: A Method for the Rapid Estimation of
% a Visual Evoked Potential" Neuroimage 32: 1549?1561.
%
% Created: 2019-03-01 Scott Bressler Apple, Inc.
X = fft(x,2^nextpow2(2*size(x,1)-1));
R = ifft(abs(X).^2);
N = length(x); % number of time samples in x
R = R./N; % Biased autocorrelation estimate
rxx = toeplitz(R(1:length(x)),conj(R(1:length(x))));

rxy = real(ifft(fft(x).*conj(fft(y))))./N;

w = flipud(rxx*rxy);