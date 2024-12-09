% statistique max global

close all;

addpath('./TF-Toolbox/TF_Toolbox');
addpath('./TF-Toolbox/TF_Toolbox/fig/');

Nfft = 512;
L = 512;
t = (0:(L - 1))/L;

B = 500;
sigma_w = 1/sqrt(B);

g = gauss_win(L, sigma_w);
g_L2Norm = sqrt(sum(g.^2));

%% realization
N_rep = 5000;
N_gumbel = L*Nfft;

l_spectro = 1/2;
c_ln = log(N_gumbel*40)/l_spectro;

% v_exp_max = [];
for r=1:N_rep
    if mod(r, 10) == 0
        fprintf("%u/%u\n", r, N_rep);
    end
    
    nri_std = 1/sqrt(2);
    n_var = 2*(nri_std^2);
    
    noise = nri_std*randn(1, L) + nri_std*1i*randn(1, L);
    STFT = stft(noise, Nfft, g);
    NSpectr = abs(STFT).^2/(nri_std^2*g_L2Norm^2);
    
    v_exp_max = [v_exp_max; max(max(NSpectr))];
end

load('data_global_max_5000_512p512.mat');

v_center1 = v_exp_max - c_ln;
g_max = max(v_center1);
g_min = min(v_center1);

% return;

%% Generate reference for comparison
%  cf. ON SPECTROGRAM LOCAL MAXIMA, P. Flandrin

% v_exp_max2 = max(exprnd(1/l_spectro, N_gumbel, N_rep));
% v_center2 = v_exp_max2 - c_ln;
% g_max = max(v_center2);
% g_min = min(v_center2);

g_sig = 1/l_spectro;
g_mu = 0;

x = g_min:(g_max - g_min)/500:g_max;

z = -(x - g_mu)/g_sig;
g_pdf = 1/g_sig*exp(z - exp(z));

%% figures

figure;
histogram(v_center1, 'Normalization', 'pdf');
hold on;
plot(x, g_pdf, 'r-');
hold off;