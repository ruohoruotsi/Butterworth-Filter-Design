%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% IOHAVOC: test drive butter unit tests cases 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% load up test wave files
% x = wavread('/Users/iorif/Desktop/mono_noise_test.wav');
x = wavread('/Users/iorif/Desktop/SinWave.wav');
x = x - mean(x);        % remove DC
x = x / max(abs(x));    % normalize

% generate Analouge lowpass prototype
[z, p, k] = buttap(8)

fp1 = 2*pi*800.2;
fp2 = 2*pi*1531.2;

% generate Analogue lowpass coefficients
% [z, p, k] = butter(8, 500/22050, 's');
% [z, p, k] = butter(8, fp1, 's')

% generate Analogue highpass coefficients
% [z, p, k] = butter(8, fp1, 'high', 's');

% generate Analogue bandpass coefficients
%[z, p, k] = butter(8, [fp1 fp2], 'bandpass', 's')

% generate Analogue bandstop coefficients
% [z, p, k] = butter(8, [fp1 fp2], 'stop', 's')


% do bilinear transform
[Zd, Pd, Kd] = bilinear(z, p, k, 48000 )
[sos,g] = zp2sos(Zd,Pd,Kd)


%%---
y = sosfilt(sos, x);
plot(real(g)*y)

figure(2);
plot(x)

%%---
% nsamps = 256;
% x = real(g)*[1,zeros(1,nsamps-1)];
% Bs = sos(:,1:3); % Section numerator polynomials
% As = sos(:,4:6); % Section denominator polynomials
% [nsec,temp] = size(sos);
% for i=1:nsec
%    x = filter(Bs(i,:),As(i,:),x); % Series sections
% end
% 
% 
% plot(x)


% plot crap
fs = 48000;
Bs = sos(:,1:3); % Section numerator polynomials
As = sos(:,4:6); % Section denominator polynomials
[nsec,temp] = size(sos);
nsamps = 256;
x = g*[1,zeros(1,nsamps-1)];
x = sosfilt(sos, x);

for i=1:nsec
  x = filter(Bs(i,:),As(i,:),x); % Series sections
end


figure(3);
X=fft(x); % sampled frequency response
f = [0:nsamps-1]*fs/nsamps; grid('on');
axis([0 fs/2 -100 5]); legend('off');
plot(f(1:nsamps/2),20*log10(X(1:nsamps/2)));




%
% fpHz = 500
% fsHz = 44100
% % prewarping tests
% fp = 2*pi*fpHz 
% fs = fpHz/tan(fpHz/fsHz/2)
% 
% % Do bilinear transformation
% pd = (1+p/fs)./(1-p/fs)                
% zd = (1+z/fs)./(1-z/fs)
% % kd = real(k*prod(fs-z)./prod(fs-p));
% 


% clear
% % load up test wave files
% x = wavread('/Users/iorif/Desktop/mono_sine_test.wav');
% x = x - mean(x);        % remove DC
% x = x / max(abs(x));    % normalize
% 
% 
% fc = 1500; % Cut-off frequency (Hz)
% fs = 44100; % Sampling rate (Hz)
% order = 8; % Filter order
% [B,A] = butter(order,2*fc/fs); % [0:pi] maps to [0:1] here
% [sos,g] = tf2sos(B,A)
% 
% % sos =
% %  1.00000  2.00080   1.00080  1.00000  -0.92223  0.28087
% %  1.00000  1.99791   0.99791  1.00000  -1.18573  0.64684
% %  1.00000  1.00129  -0.00000  1.00000  -0.42504  0.00000
% % 
% % g = 0.0029714
% %
% % Compute and display the amplitude response
% 
% Bs = sos(:,1:3); % Section numerator polynomials
% As = sos(:,4:6); % Section denominator polynomials
% [nsec,temp] = size(sos);
% nsamps = 256; % Number of impulse-response samples
% % Note use of input scale-factor g here:
% % x = g*[1,zeros(1,nsamps-1)]; % SCALED impulse signal
% 
% for i=1:nsec
%   x = filter(Bs(i,:),As(i,:),x); % Series sections
% end
% 
% %
% % figure(2);
% % X=fft(x); % sampled frequency response
% % f = [0:nsamps-1]*fs/nsamps; grid('on');
% % axis([0 fs/2 -100 5]); legend('off');
% % plot(f(1:nsamps/2),20*log10(X(1:nsamps/2)));

