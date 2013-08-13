%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% IOHAVOC: butter heq tests
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;

% -----------------------------------------------------------------------------

disp(' ');
disp('Butterworth boost example:');
disp('--------------------------');
disp('N = 4; fs=40; f0 = 5; Df = 4; w0 = 2*pi*f0/fs; Dw = 2*pi*Df/fs;'); 
disp('G0 = 0; G = 12; GB = 9; type = 0;');
disp('[B,A,Bh,Ah] = hpeq(N, G0, G, GB, w0, Dw, type);');
disp('f = linspace(0,20,1001); w = 2*pi*f/fs;');
disp('H = 20*log10(abs(fresp(B,A,w)));');
disp('[w1,w2] = bandedge(w0,Dw); f1 = fs * w1/2/pi; f2 = fs * w2/2/pi;');
disp('plot(f,H,''r-'', [f1,f2],20*log10([GB,GB]),''b.'');');



% ----------------------------------------------------------------------
% setup arguments
% ----------------------------------------------------------------------
% f1 = .8
% f2 = 3

N  = 8;             % filter order
fs = 40;            % sample rate


%Highshelf:
%f0 = 20; Df = f0 - 13;

%lowshelf
%f0 = 0; Df = .8;

% parametric
%f0 = 5; Df = 13;

%%%%%%%%%%%%%%%%%
f1 = .5
f2 = 2
f0 =(f1 + f2)/2    % center freq for parametric
Df = f2 - f1;       % delta_freq -> bandwidth


w0 = 2*pi*f0/fs; 
Dw = 2*pi*Df/fs

% gain
G0 = 0; % gain - reference
G = 12; % gain - bandwidth boost/cut
GB = G*.75; % gain - boost/cut


% ----------------------------------------------------------------------
% setup arguments
% ----------------------------------------------------------------------
% [B,A,Bh,Ah] = hpeq1(N, G0, G, GB, w0, Dw, type)

% function [B,A,Bh,Ah] = hpeq1(N,G0,G,GB,w0,Dw,type)

G0 = 10^(G0/20); 
G  = 10^(G/20); 
GB = 10^(GB/20); 


r = rem(N,2); L = (N-r)/2;

Bh = [1 0 0]; 
Ah = [1 0 0]; 

B = [1 0 0 0 0]; 
A = [1 0 0 0 0];

if G==G0, return; end              % no filtering if G=G0

c0 = cos(w0); 

if w0==0,    c0=1;  end    % special cases
if w0==pi/2, c0=0;  end
if w0==pi,   c0=-1; end

WB = tan(Dw/2);
e = sqrt((G^2 - GB^2)/(GB^2 - G0^2)); 

g = G^(1/N); g0 = G0^(1/N); 

%% Butterworth
b = WB / e^(1/N);

% if r==1,
%     D = b + 1;
%     Bh(1,:) = [g*b+g0, g*b-g0, 0]/D;
%     Ah(1,:) = [1, (b-1)/D, 0];
% 
%     B(1,:) = [g0+g*b, -2*g0*c0, g0-g*b, 0, 0]/D;
%     A(1,:) = [1, [-2*c0, 1-b, 0, 0]/D];
% end

for i=1:L,
    phi = (2*i-1)*pi/(2*N);
    si = sin(phi);
    D = b^2 + 2*b*si + 1;

    b0h = (g^2*b^2 + 2*g0*g*b*si + g0^2)/D;        % 2nd order sections
    b1h = 2*(g^2*b^2 - g0^2)/D;
    b2h = (g^2*b^2 - 2*g0*g*b*si + g0^2)/D;
    a1h = 2*(b^2 - 1)/D;
    a2h = (b^2 - 2*b*si + 1)/D;
    Bh(i+1,:) = [b0h, b1h, b2h];
    Ah(i+1,:) = [ 1,  a1h, a2h];

    b0 = (g^2*b^2 + g0^2 + 2*g*g0*si*b)/D;          % 4th order sections
    b1 = -4*c0*(g0^2 + g*g0*si*b)/D;
    b2 = 2*((1+2*c0^2)*g0^2 - g^2*b^2)/D;
    b3 = -4*c0*(g0^2 - g*g0*si*b)/D;
    b4 = (g^2*b^2 + g0^2 - 2*g*g0*si*b)/D;
    a1 = -4*c0*(1 + si*b)/D;
    a2 = 2*(1+2*c0^2 - b^2)/D;
    a3 = -4*c0*(1 - si*b)/D;
    a4 = (b^2 - 2*si*b + 1)/D;
    B(i+1,:) = [b0, b1, b2, b3, b4];
    A(i+1,:) = [1,  a1, a2, a3, a4];
end

%     Bh(:,2) = c0*Bh(:,2);	        % change sign if w0=pi
%     Ah(:,2) = c0*Ah(:,2);

% if c0==1 | c0==-1 	        % LP or HP shelving filter
%    B = Bh;                      % B,A are second-order
%    A = Ah;
%    B(:,2) = c0*B(:,2);	        % change sign if w0=pi
%    A(:,2) = c0*A(:,2);
%    B(:,4:5) = 0;                % make them (L+1)x5
%    A(:,4:5) = 0;                % for convenience in using fresp
% end  




% ----------------------------------------------------------------------
% setup frequency response
% ----------------------------------------------------------------------

f = linspace(0,20,1001); w = 2*pi*f/fs;

H = 20*log10(abs(fresp(B,A,w)));         %<---- CHANGE sh->eq!!
%H = 20*log10(abs(fresp(Bh,Ah,w,w0)));  % alternative computation of frequency response

[w1,w2] = bandedge(w0,Dw); f1 = fs * w1/2/pi; f2 = fs * w2/2/pi;




% ----------------------------------------------------------------------
% Display
% ----------------------------------------------------------------------


plot(f,H,'r-', [f1,f2],([GB,GB]),'b.');
ylim([-8 14]); ytick(-6:3:12);
xlim([0,20]); xtick(0:2:20);
title('Butterworth Boost {\it N} = 4');
xlabel('{\it f}  (kHz)'); ylabel('dB');
grid;



% ----------------------------------------------------------------------
% first convert equalizer fourth order sections
% ----------------------------------------------------------------------
Bh
Ah

B
A

[SOSeq, Geq] = tf2sos(B,A);


% ----------------------------------------------------------------------
% filter actual sound
% ----------------------------------------------------------------------


% plot crap
fs = 40000;
Bs = SOSeq(:,1:3); % Section numerator polynomials
As = SOSeq(:,4:6); % Section denominator polynomials
[nsec,temp] = size(SOSeq);
nsamps = 256;
x = Geq*[1,zeros(1,nsamps-1)];


% Filter using direct form transposed II
for i=1:nsec/2
  %x = filter(B(i,:),A(i,:),x); % Series sections %<---- CHANGE sh->eq!!
  x = filter(Bh(i,:),Ah(i,:),x); % WORKS FOR SHELVING!!!
end


figure(2);
X=fft(x); % sampled frequency response
f = [0:nsamps-1]*fs/nsamps; grid('on');
axis([0 fs/2 -100 5]); legend('off');
plot(f(1:nsamps/2),20*log10(X(1:nsamps/2)));

