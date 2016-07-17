clear all;
close all;

Tx=1; % period of input signal in sec
W=1024; % Nyquist rate (Hz); W/2 Hz is the signal bandlimit
oversample=1; % oversampling factor (wrt the Nyquist rate) of the simulated input signal
N=Tx*W*oversample; % length of simulated input signal (when oversample=1 there are N Nyquist periods in Tx)
t=0:1/(W*oversample):Tx-1/(W*oversample); % observation interval in sec [0,Tx]
K=1:1:100; % sparsity parameter; number of frequencies present in input signal
sig = [];
for i=1:1:100
    for ii = 1:30
       [x,X,freq]=tonesparse(t,Tx,N,K(i)); % generate spectrally-sparse multitone signal 
       sig = [sig x];
    end
end
%% 使用KSVD算法进行字典训练并写入xlsx文件
codebook_size = 2000;
errGoal = 1e-6;
[Dic,co]= K_SVD(sig,codebook_size,errGoal);
xlswrite('Dic.xlsx',Dic);