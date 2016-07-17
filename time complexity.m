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


Mt=[32,64,128,256,512,1024];
t1=zeros(length(Mt),1);t2=t1;t3=t1;
for i=1:length(Mt)
    M=Mt(i);

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



 
%for cnt = 1:10
       % SNR = SNRmtx;
        
        tic;
        [x_hat1,X_hat1,freq_hat1]=rd_recovery1(y3,t,S,N,M,Tx,K); % recover input signal rd_recovery4(y,t,S,N,M,Tx)
        t1(i)=toc;
       
        tic;
        [x_hat2,X_hat2,freq_hat2]=rd_recovery4(y3,t,S,N,M,Tx); % recover input signal rd_recovery4(y,t,S,N,M,Tx)
        t2(i)=toc;
        
        tic;
        [x_hat3,X_hat3,freq_hat3]=rd_recovery8(y3,t,S,N,M,Tx); % recover input signal rd_recovery4(y,t,S,N,M,Tx)
        t3(i)=toc;
          
        
       

end



plot(Mt/N,t1,'r-o','Linewidth',2) %plot(x,y,'Linewidth',3)
hold on
plot(Mt/N,t2,'b-^','Linewidth',2) %plot(x,y,'Linewidth',3)
hold on
plot(Mt/N,t3,'k-s','Linewidth',2) %plot(x,y,'Linewidth',3)
grid on
legend('OMP','SAMP','improved-SAMP')
xlabel('compression rate'),ylabel('time complexity')
%title('RD系统算法运行时间复杂度随压缩率变化情况')



