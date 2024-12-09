close all

addpath('./toolbox/TF_Toolbox');
addpath('./toolbox/TF_Toolbox/fig/');

%bat signal  
load -ascii batsig.txt
s_origin = batsig;
s = s_origin(145:end);
s = s(:).';

s = hilbert(s);

sigma_w = 0.13;

s_ref = s;
IFs = [0, 512];
Amp = ones(size(s));
s_name = 'bat';
NR = 3;

s_all = s;
Lx = length(s_all);
Fs = 160000;
Nfft = 1024;

Tx = (0:(Lx - 1))/Fs;
Fx = (0:(Nfft - 1))*Fs/Nfft;
% [STFT, T] = sstn(s_all, sigma_w, Nfft);
[g, ~] = gauss_win(Lx, sigma_w);
V = stft(s_all, Nfft, g);

% Tx = (0:(Lx-1))/Lx;
% F_vec = (0:N_f-1)*Lx/Nfft;

% xi0 = min(IFs) - 100;
% xi1 = max(IFs) + 100;
% F2 = F_vec.*(F_vec > xi0).*(F_vec < xi1);

TFRsc(Tx*1000, Fx/1000, abs(V), 'xunit', 'ms', 'yunit', 'kHz');
ylim([0, 78]);
write_figfiles('bat');
