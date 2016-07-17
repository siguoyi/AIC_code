 %% Sampling spectrally-sparse multitone signals with the Random Demodulator

% This top-level script (m-file) demonstrates the sampling and recovery of sparse
% multitone signals using the Random Demodulator (RD). The RD was originally proposed by
%  J.A. Tropp, J.N. Laska, M.F. Duarte, J.K. Romberg, and R.G. Baraniuk in "Beyond Nyquist:
%  Efficient Sampling of Sparse Bandlimited Signals", Information Theory, IEEE Trans. on,
%  vol. 56, no. 1, pp. 520-544, Jan 2010.

% When the script is executed MATLAB will generate a discrete multitone signal 
% (simulating a continuous one), sample it, and recover the original signal from 
% its samples. The modeling and system parameters are set in this script. For a 
% complete description of the parameters, refer to the technical report,
% "Sampling Sparse Multitone Signals with a Random Demodulator", accompanying 
% the CTSS Sampling Toolbox or the above reference.

% The sub-routine rd_recovery.m uses the software package Sparsify 0.4 by Thomas
% Blumensath. Available at: www.personal.soton.ac.uk/tb1m08/sparsify/.

% Copyright (c) 2010 by Michael A Lexa, Mike E Davies, John S Thompson
% Institute of Digital Communications, University of Edinburgh
%
% This work is licenced under the Creative Commons Attribution 3.0
% Unported License. To view a copy of this licence, visit
% http://creativecommons.org/licenses/by/3.0/ or send a letter to
% Creative Commons, 171 Second Street, Suite 300, San Francisco,
% California 94105, USA.

clear
addpath(genpath('sparsify_0_4'));

%% Parameters for input signal
Tx=1; % period of input signal in sec
W=1024; % Nyquist rate (Hz); W/2 Hz is the signal bandlimit
oversample=1; % oversampling factor (wrt the Nyquist rate) of the simulated input signal
K=10; % sparsity parameter; number of frequencies present in input signal
N=Tx*W*oversample; % length of simulated input signal (when oversample=1 there are N Nyquist periods in Tx)
t=0:1/(W*oversample):Tx-1/(W*oversample); % observation interval in sec [0,Tx]

%% Parameters for random demodulator
% M = 64*[1:2:8];
m=2;
M=64; % the ratio M/N specifies rate of sampling process relative to Nyquist rate 
L=N; % there are L Nyquist periods in one period of p(t) (for Tropp's analysis Tx=Tp or L=N)
h=[ones(1,N/M), zeros(1,(M-1)*N/M)]; % impulse response for ideal integrator
hf = abs(fft(h))*2/(N-1);
[x,X,freq]=tonesparse(t,Tx,N,K); % generate spectrally-sparse multitone signal
%x=x+sqrt(.01)*randn(size(x)); % uncomment to add white Gaussian signal noise

S=2*randi([0,1],L,m)-1; % generate random +/- 1 sequences
%S=pnmat(L,1); % generate random +/- 1 sequences
y1=zeros(N,m);
y3=zeros(M/m,m);
for channel=1:m
    y1(:,channel)=x.*S(:,channel);
    y2(:,channel)=conv(h,y1(:,channel));
    y4=y2(:,channel);
    y3(:,channel)=downsample(y4((N/M)*m:N),(N/M)*m);
end
%y1=x.*S; % multiply input signal by random +/- 1 sequence 
yy2=conv(h,y1(:,unidrnd(m))); % filter demodulated signal with ideal integrator
yy3=downsample(yy2((N/M)*m:N),(N/M)*m); % sample
yf1=abs(fft(y1(:,unidrnd(m))))*2/N;
yf2=yf1'.*hf; 
figure(1)
subplot(4,1,1)
plot(t,x)
title('多音信号时域波形')
subplot(4,1,2)
plot(t,y1(:,unidrnd(m)))
title('调制信号时域波形')
subplot(4,1,3)
plot(t,yy2(1:N))
title('滤波信号时域波形')

%[y,S]=rd_sampling(x,L,h,N,M); % random demodulation sampling
[x_hat,X_hat,freq_hat]=rd_recovery6(y3,t,S,N,M,Tx,K,m); % recover input signal
subplot(4,1,4)
plot(t,x_hat)
title('重构信号时域波形')
figure(5)
subplot(4,1,1)
ft=(0:length(t)-1)*W/length(t)-W/2;
plot(ft,abs(fftshift(abs(fft(x))*2/N)));
title('多音信号频域波形');
subplot(4,1,2)
plot(ft,abs(fftshift(abs(fft(yf1))*2/N)));
title('调制信号频域波形');
subplot(4,1,3)
plot(ft,abs(fftshift(yf2)));
title('滤波后信号频域波形');
subplot(4,1,4)
plot(ft,abs(fftshift(abs(fft(x_hat))*2/N)));
title('重构信号频域波形');

%% Results
fprintf('Nyquist rate %3.3f Hz\n',W);
fprintf('System sampling rate %3.3f Hz\n',(M*W)/N);
fprintf('Subsampling factor %3.3f\n',M/N);
fprintf('Original Frequencies (Hz)\n');
fprintf('%6.10f\n',sort(-freq/Tx));
fprintf('Recovered Frequencies (Hz)\n');
fprintf('%6.2f\n',freq_hat/Tx);
fprintf('Original FS Coefficients (real)\n');
fprintf('%6.2f\n',flipud(X(X>0)));
fprintf('Recovered FS Coefficients (real)\n');
fprintf('%6.2f\n',X_hat(X_hat>0));
%fprintf('Average squared error =%2.10f\n',sum((1/length(x))*(real(x)-real(x_hat)).^2));
fprintf('Average squared error =%2.10f\n',20*log10(norm(x,2)/norm((x-x_hat),2)));
%fprintf('Average squared error2 =%2.10f\n',20*log10(norm(real(x),2)/norm((real(x)-real(x_hat)),2)));
%osnr1=20*log10(norm(x,2)/norm((x-r'),2));
figure(2) % graphically display recovered frequencies and FS coefficients compared to true values
stem(-N/2:N/2-1,circshift(flipud(X),1)); hold on; 
stem(-N/2:N/2-1,real(X_hat),'r','Marker','*');
xlabel('Frequency Index'); ylabel('Magnitude'); title('FS Coefficients (Magnitude)'); 
legend('Original FS coefficients','Recovered FS coefficients')

figure(3) % plot fft of original multitone signal x and reconstructed signal x_hat
f= (-length(x)/2:length(x)/2-1)./length(x) * (W*oversample);
h1=subplot(211); stem(f,abs(fftshift(abs(fft(x))*2/N))); 
title('Fourier Spectrum--Input signal'); xlabel('Hz'); ylabel('Magnitude')
h2=subplot(212); stem(f,abs(fftshift(abs(fft(x_hat))*2/N))); 
title('Fourier Spectrum--Reconstructed signal');xlabel('Hz'); ylabel('Magnitude')
linkaxes([h1 h2])

figure(4) % plot original multitone signal x and reconstructed signal x_hat
h3=subplot(311); stem(t,real(x),'Marker','none'); 
xlabel('Seconds'); title('Original Signal (real part)');
h4=subplot(312); stem(t,real(x_hat),'r','Marker','none'); 
xlabel('Seconds'); title('Reconstructed Signal (real part)');
h5=subplot(313); stem(t,real(x)-real(x_hat),'Marker','none'); 
xlabel('Seconds'); title('Difference Signal');
linkaxes([h3 h4 h5])