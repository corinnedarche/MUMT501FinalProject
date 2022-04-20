% Script name: onsetImplementation.m
% Created by: Corinne Darche
% Description: This script is an implementation of Depalle and Scherrer's
% 2015 paper "Onset Time Estimation for the Exponentially Damped Sinusoids
% Analysis of Percussive Sounds"

% This script was made as part of a final project on Onset Detection for
% MUMT501: Digital Signal Processing (McGill University, Winter 2022)

% Define variables using the ones mentioned in the paper
clear

Fs = 44100;
T = 1/Fs;
N = 1024; % FFT size
H = 512; % Hop size
gamma_rough = 0.3;
tau_rough = 0.15; % absolute threshold
p = 5; % order of median filter
l = 0.5; % how much is absolute threshold affected by d_(median,p)
alpha = 6; % threshold, in dB
I_fine = 900; % time interval
I_rough = 2205;

%% Create the synthetic signal

% Synthetic sound: exponentially damped sinusoid
t = 0:T:1;

p1 = 0.7.*exp(-t.*30).*cos(2*pi*261.63.*t + pi) + 0.68.*exp(-t.*30.4).*cos(2*pi*263.98.*t + 0.95*pi);
p2 = 0.4.*exp(-t.*45).*cos(2*pi*261.63.*t + pi/2) + 0.49.*exp(-t.*45.9).*cos(2*pi*263.98.*t + 0.789*(pi/2));
p3 = 0.3.*exp(-t.*60).*cos(2*pi*261.63.*t + pi) + 0.35.*exp(-t.*60.2).*cos(2*pi*263.98.*t + 0.95*pi);
p4 = 0.26.*exp(-t.*75).*cos(2*pi*261.63.*t + pi/2) + 0.269.*exp(-t.*75.5).*cos(2*pi*263.98.*t + 0.789*(pi/2));
p5 = 0.25.*exp(-t.*90).*cos(2*pi*261.63.*t + pi) + 0.2507.*exp(-t.*90.7).*cos(2*pi*263.98.*t + 0.95*pi);

event = p1+p2+p3+p4+p5;
tran1 = [zeros(1,1000) event(1:Fs-999)];
tran2 = [zeros(1,Fs/2) event(1:(Fs/2)+1)];

x = tran1 + tran2;

%% Read in an audio file (taken from Leveau dataset)

[x, Fs] = audioread("Leveau/sounds/guitar2.wav");
x = x.';
T = 1/Fs;
info = audioinfo("Leveau/sounds/guitar2.wav");
dur = info.Duration;
t = 0:T:dur;
t(end) = [];
load("Leveau/goodlabels/guitar2.mat")

%% First onset determination: STFT with "rough" time resolution

[X, fSTFT, tSTFT] = stftAlt(x,N,H,N,Fs); % Use STFT implementation that takes N and H as input

% Frequency-domain detection function d_f[l]
dFreq = zeros(1,size(X,2)-1);
for k=2:(size(dFreq,2))
    total = 0;
    for m=1:(N/2)
        val = (abs(X(m,k)) - abs(X(m,k-1)))^2;
        total = total + val;
    end
    dFreq(k) = sqrt(total);
end

% Post-processing and peak-detection
dFreqNorm = postprocess(dFreq,gamma_rough);
[peaksFreq, indFreq, freqTime] = parabolicFreq(dFreqNorm,alpha,tSTFT);

threshAdaptFreq = tau_rough + l.*medfilt1(dFreqNorm,p);

onsetRoughTime = freqTime;
onsetRoughPeaks = peaksFreq;

for i=1:size(indFreq)
    if peaksFreq(i) < threshAdaptFreq(indFreq(i))
        onsetRoughTime(i) = [];
        onsetRoughPeaks(i) = [];
    end
end

ITimeRough = I_rough * T;

if size(onsetRoughTime,2) > 1
    k = 1;
    while k < size(onsetRoughTime, 2) - 1
        while abs(onsetRoughTime(1,k) - onsetRoughTime(1,k+1)) < ITimeRough && k < (size(onsetRoughTime, 2)-1)
            if onsetRoughPeaks(k) > onsetRoughPeaks(k+1)
                onsetRoughTime(k+1) = [];
            else
                onsetRoughTime(k) = [];
            end
        end
        k = k + 1;
    end
end
%% Second onset determination: finer time resolution

% Define additional variables for finer resolution

J = 400;
v = 10^(-4);
gamma_fine = 0.1;
tau_fine = 0.5;

nMax = size(x,2) - J;

% Time-domain Function

dTime = zeros(1,(nMax - J)+1);
a = 1;

for ind=(J+1):nMax
    num = 0;
    denom = 0;
    additional = 0;
    for m=(ind+1):(ind+J)
        num = num + (x(m))^2;
    end
    for o=(ind-J+1):(ind-1)
        denom = denom + ((x(o))^2 + v);
    end
    for k=(ind+1):(ind+J)
        additional = additional + (x(k))^2;
    end
    dTime(a) = (1/J)*log10(num/(denom)) * additional;
    a = a+1;
end

dTime(isnan(dTime)) = 0;

tTime = (J:nMax).*T;

dTimeNorm = postprocess(dTime,gamma_fine);

threshAdaptTime = tau_fine + l.*medfilt1(dTimeNorm,p);

[peaksTime, indTime, timeTime] = parabolicTime(dTimeNorm, alpha,tTime);

onsetFineTime = timeTime;
onsetFinePeaks = peaksTime;

for i=1:size(indTime)
    if peaksTime(i) < threshAdaptTime(indTime(i))
        onsetFineTime(i) = [];
        onsetFinePeaks(i) = [];
    end
end

ITimeFine = I_fine * T;

if size(onsetFineTime,2) > 1
    k = 1;
    while k < (size(onsetFineTime, 2) - 1)
        while abs(onsetFineTime(k) - onsetFineTime(k+1)) < ITimeFine && k < (size(onsetFineTime, 2)-1)
            if onsetFinePeaks(k) > onsetFinePeaks(k+1)
                onsetFineTime(k+1) = [];
            else 
                onsetFineTime(k) = [];
            end
        end
        k = k + 1;
    end
end

%% Plotting onsets

figure(1)

plot(t,x)
for k = 1:size(labels_time, 2) xline(labels_time(1,k), 'g-'); end

%legend('Synthetic Signal', 'Rough Onset', 'Refined Onset')

xlabel('Time (s)')
ylabel('Amplitude')



