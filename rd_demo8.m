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
W=2048; % Nyquist rate (Hz); W/2 Hz is the signal bandlimit
oversample=1; % oversampling factor (wrt the Nyquist rate) of the simulated input signal
K=10; % sparsity parameter; number of frequencies present in input signal
N=Tx*W*oversample; % length of simulated input signal (when oversample=1 there are N Nyquist periods in Tx)
t=0:1/(W*oversample):Tx-1/(W*oversample); % observation interval in sec [0,Tx]


CR=[2,4,8,16,32,64,128,256];
successrate=zeros(length(CR),1);
for i=1:length(CR)
    M=N/CR(i);

%% Parameters for random demodulator
% M = 64*[1:2:8];
%M=64; % the ratio M/N specifies rate of sampling process relative to Nyquist rate 
L=N; % there are L Nyquist periods in one period of p(t) (for Tropp's analysis Tx=Tp or L=N)
h=[ones(1,N/M), zeros(1,(M-1)*N/M)]; % impulse response for ideal integrator
hf = abs(fft(h))*2/(N-1);
[x,X,freq]=tonesparse(t,Tx,N,K); % generate spectrally-sparse multitone signal
%x=x+sqrt(.01)*randn(size(x)); % uncomment to add white Gaussian signal noise

S=2*randi([0,1],L,1)-1; % generate random +/- 1 sequences
%S=pnmat(L,1); % generate random +/- 1 sequences
y1=x.*S; % multiply input signal by random +/- 1 sequence 
y2=conv(h,y1); % filter demodulated signal with ideal integrator
y3=downsample(y2(N/M:N),N/M); % sample
% yf1=abs(fft(y1))*2/N;
% yf2=yf1'.*hf; 
% figure(1)
% subplot(4,1,1)
% plot(t,x)
% title('多音信号时域波形')
% subplot(4,1,2)
% plot(t,y1)
% title('调制信号时域波形')
% subplot(4,1,3)
% plot(t,y2(1:N))
% title('滤波信号时域波形')

%[y,S]=rd_sampling(x,L,h,N,M); % random demodulation sampling


success=zeros(5,1);
for cnt = 1:length(success)
       % SNR = SNRmtx;
        
        
        [x_hat1,X_hat1,freq_hat1]=rd_recovery1(y3,t,S,N,M,Tx,K); % recover input signal rd_recovery4(y,t,S,N,M,Tx)
        if length(freq_hat1)==10
            freq_hat=freq_hat1;
        else
            freq_hat=zeros(10,1);
            freq_hat(1:length(freq_hat1))=freq_hat1;
        end
            
        if freq_hat==sort(-freq)
            success(cnt,1)=1;
        else
            success(cnt,1)=0;
        end
        
       
end
   
successrate(i,1)= 100*sum(success)/(length(success));

end
% subplot(4,1,4)
% plot(t,x_hat)
% title('重构信号时域波形')
% figure(5)
% subplot(4,1,1)
% ft=(0:length(t)-1)*W/length(t)-W/2;
% plot(ft,abs(fftshift(abs(fft(x))*2/N)));
% title('多音信号频域波形');
% subplot(4,1,2)
% plot(ft,abs(fftshift(abs(fft(y1))*2/N)));
% title('调制信号频域波形');
% subplot(4,1,3)
% plot(ft,abs(fftshift(yf2)));
% title('滤波后信号频域波形');
% subplot(4,1,4)
% plot(ft,abs(fftshift(abs(fft(x_hat))*2/N)));
% title('重构信号频域波形');

plot(CR,successrate,'r-o','Linewidth',1.5) %plot(x,y,'Linewidth',3)
grid on
legend('OMP')
xlabel('compression rate'),ylabel('success percentage')
title('RD系统支撑集重构成功率随压缩率变化情况')


