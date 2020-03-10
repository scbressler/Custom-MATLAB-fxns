%--------------------------------------------------------------------------
% [yt,t] = cosineToneComplex(fs,dur,freqs,phases,amp)
%--------------------------------------------------------------------------
% Creates a tone complex in cosine-phase.
%
% INPUT VARIABLES
%   fs = sampling frequency
%   dur = duration of tone complex in seconds
%   freqs = vector of frequency components
%   phases = vector of starting phases in radians
%   amp = vector of component amplitudes
%
% OUTPUT VARIABLES
%   yt = sine tone complex vector
%    t = time/sample in seconds
%--------------------------------------------------------------------------
% Created by: Scott Bressler
% 9 AUGUST 2011
% Modified:
% 22 AUGUST 2014, SCB: added 'amp' option
% 28 OCTOBER 2014, SCB: 'amp' can have different values per harmonic
% 16 APRIL 2015, SCB: output, yt, normalized based on number of harmonics
%--------------------------------------------------------------------------
function [yt,t] = cosineToneComplex(fs,dur,freqs,phases,amp)

if nargin<4
    phases = zeros(size(freqs));
    amp = ones(size(freqs));
    disp(sprintf('Starting phases of all components = 0'));
    disp(sprintf('Amplitude of all components = 1'));
end

if(length(phases)==1)
    phases = repmat(phases,size(freqs));
end

if(length(amp)==1)
    amp = repmat(amp,size(freqs));
end

t = [0:1/fs:dur]';

phases = repmat(phases(:)',length(t),1);
amp = repmat(amp(:)',length(t),1);

cosinArg = 2*pi*freqs(:)';

% xt = cos(t*cosinArg+phases).*repmat(amp,length(t),1);
xt = amp.*cos(t*cosinArg+phases);

yt = sum(xt,2)./length(freqs);


