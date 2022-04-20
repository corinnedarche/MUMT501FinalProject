function [peakVals, ind, timeInd] = parabolicFreq(input, a, t)
% Function: parabolicFreq.m
% Created by: Corinne Darche
% Use parabolic interpolation to find the peaks with frequency-domain
% detection function as input

% Taken from MUMT307 Parabolic Interpolation code
ind = [];

threshold = 10^(a/20); % Linear gain corresponding to alpha dB (6)

% Get the indices for the peaks
for n = 2:length(input) - 2
      if input(n+1) <= input(n) && input(n-1) < input(n) && input(n) > threshold
          ind = [ind, n];
      end
end

timeInd = zeros(size(ind));

% Now use parabolic interpolation
peakVals = zeros( size( ind ) );

for n = 1:length(ind)
    alpha = input(ind(n)-1);
    beta = input(ind(n));
    gamma = input(ind(n)+1);
    Xiind = 0.5*(alpha-gamma) / (alpha - 2*beta + gamma);
    peakVals(n) = beta - 0.25 * (alpha-gamma) * Xiind;
    if Xiind > 0
        timeInd(n) = t(ind(n)) + Xiind*(t(ind(n)+1)-t(ind(n)));
    else
        timeInd(n) = t(ind(n)) + Xiind*(t(ind(n))-t(ind(n)-1));
    end
end


end

