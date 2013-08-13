%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% IOHAVOC: coeff tests for biquad BPF w/ constant skirt
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boost only filter

% fc = .2;
% fs = 1;
% bw = fs/20;
% gain = 2.0;
% 
% 
% Q   = fs/bw;
% wcT = 2*pi*fc/fs;
% 
% K = tan(wcT/2);
% V = gain;
% 
% b0 =  1 + V*K/Q + K^2;
% b1 =  2*(K^2 - 1);
% b2 =  1 - V*K/Q + K^2;
% 
% a0 =  1 + K/Q + K^2;
% a1 =  2*(K^2 - 1);
% a2 =  1 - K/Q + K^2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LPF:        H(s) = 1 / (s^2 + s/Q + 1)
%

% Fs = 44100;
% f0 = 500;
% Q = 1/2;
% 
% w0 = 2*pi*f0/Fs;
% alpha = sin(w0)/(2*Q);
% 
% b0 =  (1 - cos(w0))/2;
% b1 =   1 - cos(w0);
% b2 =  (1 - cos(w0))/2;
% 
% a0 =   1 + alpha;
% a1 =  -2*cos(w0);
% a2 =   1 - alpha;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HPF:        H(s) = s^2 / (s^2 + s/Q + 1)
%
 
% Fs = 44100;
% f0 = 5000;
% Q = 1/2;
% 
% w0 = 2*pi*f0/Fs;
% alpha = sin(w0)/(2*Q);
% 
% b0 =  (1 + cos(w0))/2;
% b1 = -(1 + cos(w0));
% b2 =  (1 + cos(w0))/2;
% a0 =   1 + alpha;
% a1 =  -2*cos(w0);
% a2 =   1 - alpha;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BPF: H(s) = (s/Q) / (s^2 + s/Q + 1)      (constant 0 dB peak gain)
% 
% Fs = 44100;
% f0 = 15000;
% Q = 12;
% 
% w0 = 2*pi*f0/Fs;
% alpha = sin(w0)/(2*Q);
% 
% b0 =   alpha;
% b1 =   0;
% b2 =  -alpha;
% a0 =   1 + alpha;
% a1 =  -2*cos(w0);
% a2 =   1 - alpha;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BPF:  H(s) = s / (s^2 + s/Q + 1)  (constant skirt gain, peak gain = Q)
%
% Fs = 44100;
% f0 = 15000;
% Q = 12;
% 
% w0 = 2*pi*f0/Fs;
% alpha = sin(w0)/(2*Q);
% 
% b0 =   Q*alpha;
% b1 =   0;
% b2 =  -Q*alpha;
% a0 =   1 + alpha;
% a1 =  -2*cos(w0);
% a2 =   1 - alpha;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peakingEQ:  H(s) = (s^2 + s*(A/Q) + 1) / (s^2 + s/(A*Q) + 1)
% good for both boosting and peaking targeted areas, not so good match Q
% to bandwidth requirements
% 
Fs = 44100;
fl = 300;
fu = 1500;
f0 = (fu + fl)/2;
BW = log2(fu/fl)

w0 = 2*pi*f0/Fs;


% How to Q??

% Q = f0/BW
% Q = 1/ (2*sinh(log(2)/2*BW*w0/sin(w0)))
Q = sqrt(2^BW)/(2^BW - 1)
% Q = 1/(2*sinh(log(2)/2*BW)) 


dBgain = 4;

A  = 10^(dBgain/40);   %(for peaking and shelving EQ filters only)
alpha = sin(w0)/(2*Q) 

b0 =   1 + alpha*A;
b1 =  -2*cos(w0);
b2 =   1 - alpha*A;

a0 =   1 + alpha/A;
a1 =  -2*cos(w0);
a2 =   1 - alpha/A;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All Pass Filter :  H(s) = (s^2 - s/Q + 1) / (s^2 + s/Q + 1)
%  

% alpha = 1;
% 
% b0 =   1 - alpha;
% b1 =  -2*cos(w0);
% b2 =   1 + alpha;
% 
% a0 =   1 + alpha;
% a1 =  -2*cos(w0);
% a2 =   1 - alpha;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  scale coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A  = [a0 a1 a2] / a0;
B  = [b0 b1 b2] / a0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);
freqz(B,A);
title('Adjust Boost/Cut Frequency Response')
