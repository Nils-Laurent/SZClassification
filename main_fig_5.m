%% BAT SIGNAL

close all;

addpath('./TF-Toolbox/TF_Toolbox');
addpath('./TF-Toolbox/TF_Toolbox/fig/');
addpath('./project_functions');
%%

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

s_all = s;
Lx = length(s_all);
Fs = 160000;
Nfft = 1024;

c_snr = 10;
[s_noise, ~] = add_noise(s_all, c_snr);
% s_noise = s_all;

g = gauss_win(Lx, sigma_w);
STFT = stft(s_noise, Nfft, g);

[id_zeros,nv,kv] = vor_ID_zeros_bat(s_noise, sigma_w, Nfft);

id_sig = id_zeros == 1;
id_int = id_zeros == 2;

Tx=1:Lx;
Fx=1:Nfft;

TFRsc(Tx, Fx, abs(STFT));
hold on;
voronoi(nv, kv);
plot(nv(id_sig), kv(id_sig), 'rx');
plot(nv(id_int), kv(id_int), 'co');
hold off;
ylim([0, 510]);
set(gca,'XTick',[], 'YTick', []);
write_figfiles("FIG5_bat");