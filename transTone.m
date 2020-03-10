function ttone = transTone(tfs,env,fs,dur,depth)
%--------------------------------------------------------------------------
%ttone = transTone(tfs,env,fs,dur,depth)
%--------------------------------------------------------------------------
%INPUTS:
%   tfs = frequency of carrier in Hz
%   env = frequency of envelope in Hz
%    fs = sampling frequency in Hz
%   dur = tone duration in seconds
% depth = modulation depth in dB
%           0 dB = fully modulated
%          -6 dB = 50% modulated
%         -12 dB = 25% modulated
%
%OUTPUT
% ttone = transposed tone waveform
%--------------------------------------------------------------------------
% Created by: Scott Bressler
% 17-October-2013
%
% From:
%   Bernstein, LR and Trahiotis, C (2002) "Enhancing sensitivity to
%   interaural delays at high frequency by using "transposed stimuli," JASA
%   112(3), Part 1, pp 1026-1036
%
% Modified: 8-Jun-2016 to include noise carrier [SB]
%--------------------------------------------------------------------------

t = 0:1/fs:dur;     %time vector
N = length(t);      %vector length in samples
f = (0:N-1)/(N/fs); %frequency vector

%% CREATE COMPONENT WAVEFORMS
% Carrier waveform in sine phase
if(tfs==0) % TFS is BPF white noise (3-10kHz)
    tfst = randn(1,length(t));
    tfst = tfst/max(abs(tfst));
    
    Ny = fs/2;
    Fstop1 = 2500/Ny;  % First Stopband Frequency
    Fpass1 = 3000/Ny;  % First Passband Frequency
    Fpass2 = 10000/Ny;  % Second Passband Frequency
    Fstop2 = 10500/Ny;  % Second Stopband Frequency
    Astop1 = 60;      % First Stopband Attenuation (dB)
    Apass  = 1;       % Passband Ripple (dB)
    Astop2 = 60;      % Second Stopband Attenuation (dB)

h = fdesign.bandpass('fst1,fp1,fp2,fst2,ast1,ap,ast2', Fstop1, Fpass1, ...
    Fpass2, Fstop2, Astop1, Apass, Astop2);

Hd = design(h, 'equiripple');

tfst = filter(Hd,tfst);
else
    tfst = sin(2*pi*tfs*t);
end

% Envelope waveform
% Half-wave rectified
envt = max((10^(depth/20))*sin(2*pi*env*t),zeros(size(t)));
envXw = fft(envt);

cutoff = 2000;      %envelope cutoff frequency (default = 2000 Hz)

% Zero-out frequencies greater than cutoff frequency in the frequency
% domain, and inverse Fourier transform back to time domain
k1 = find(f<cutoff,1,'last');
k2 = find(envXw==conj(envXw(k1)));

envXw(1,[k1:k2]) = 0;
envt = ifft(envXw,'symmetric')+(1-(10^(depth/20)));

% Combine Carrier and Envelope and scale maximum amplitude to 0.99
ttone = (tfst.*envt)';
ttone = 0.99*(ttone/max(abs(ttone)));

