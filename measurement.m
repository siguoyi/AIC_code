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
Mtx=[32,64,128,512]; % the ratio M/N specifies rate of sampling process relative to Nyquist rate 
L=N; % there are L Nyquist periods in one period of p(t) (for Tropp's analysis Tx=Tp or L=N)

[x,X,freq]=tonesparse(t,Tx,N,K); % generate spectrally-sparse multitone signal
%x=x+sqrt(.01)*randn(size(x)); % uncomment to add white Gaussian signal noise
%noise = random('norm', 0, 0.1, [N 1]);%º”‘Î…˘
sigma=0.1;
noise = sigma*randn([N 1]);%º”‘Î…˘
%figure(6)
%plot(t,noise);
SNRmtx=-10:5:80;
NoiseEnergy = norm(noise,2)^2;
SignalEnergy = norm(x,2)^2;
CurrentSNR = 10*log10(SignalEnergy/NoiseEnergy);
CurrentSNR_val=10^(CurrentSNR/10);

%  Maximal outside iteration number:
% opts.maxit=10;
% %  Send the true solution merely for information display
% opts.xs = X;
% opts.sigma = sqrt(NoiseEnergy);

S=2*randi([0,1],L,1)-1; % generate random +/- 1 sequences

y1=x.*S; % multiply input signal by random +/- 1 sequence 
y11=noise.*S;
meanoutputSNR1=zeros(length(Mtx),length(SNRmtx));%meanoutputSNR2=meanoutputSNR1;   meanoutputSNR3=meanoutputSNR1;% meanoutputSNR3=meanoutputSNR2;
for c=1:length(Mtx)
    M=Mtx(c);


h=[ones(1,N/M), zeros(1,(M-1)*N/M)]; % impulse response for ideal integrator
hf = abs(fft(h))*2/(N-1);
y2=conv(h,y1); % filter demodulated signal with ideal integrator
y22=conv(h,y11);
y3=downsample(y2(N/M:N),N/M); % sample
y33=downsample(y22(N/M:N),N/M);






for Ki = 1:length(SNRmtx)
    
        SNR = SNRmtx(Ki);
        SNR_val = 10^(SNR/10);
        y4=y3+y33*sqrt(CurrentSNR_val/SNR_val);
 
        [x_hat1,X_hat1,freq_hat1]=rd_recovery8(y4,t,S,N,M,Tx);
        meanoutputSNR1(c,Ki) =20*log10(norm(x,2)/norm((x-x_hat1),2));
        
end
end




figure(1)

 plot(SNRmtx,meanoutputSNR1(1,:),'r-o','Linewidth',1.5);
 hold on
 plot(SNRmtx,meanoutputSNR1(2,:),'-gs','Linewidth',1.5);
 hold on
 plot(SNRmtx,meanoutputSNR1(3,:),'-b+','Linewidth',1.5);
 hold on
 plot(SNRmtx,meanoutputSNR1(4,:),'-kd','Linewidth',1.5);
 grid on
xlabel('input SNR(dB)');
ylabel('mean output SNR(dB)');
legend('M=32','M=64','M=128','M=512');
%title('SAMP Algorithm');