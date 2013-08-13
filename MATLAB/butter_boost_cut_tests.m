%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% IOHAVOC: tests for boost/cut
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% load up test wave files
% x = wavread('/Users/iorif/Desktop/mono_noise_test.wav');
x = wavread('/Users/iorif/Desktop/SinWave.wav');
x = x - mean(x);        % remove DC
x = x / max(abs(x));    % normalize


fs = 48000;
fp1 = 8000 / fs
fp2 = 15000 /fs
fc = 6500 / fs;

% generate Analogue lowpass coefficients
% [z, p, k] = butter(8, 500/22050, 's');
% [z, p, k] = butter(8, fp1, 's')

% generate Analogue highpass coefficients
[B, A] = butter(8, [fp1 fp2], 'bandpass');

% generate Analogue bandpass coefficients
%[B, A] = butter(8, [fp1 fp2], 'stop')

% generate Analogue bandstop coefficients
% [B, A] = butter(8, fc, 'low');


y = filter(B, A, x);

z = x + y;

spectrogram(z,2048, 256, 2048, 'yaxis'); 
title('Linear Chirp: start at 100Hz and cross 200Hz at t=1sec');
