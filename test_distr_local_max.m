% statistique max local

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
c_ln = 9;

vec = [];
for n_rep = 1:NRep
    if mod(r, 10) == 0
        fprintf("%u/%u\n", r, N_rep);
    end
    
    noise = randn(1, L) + 1i*randn(1, L);
    STFT = stft(noise, Nfft, g);
    NSpectr = abs(STFT).^2/(g_L2Norm^2);
    
    BW = imregionalmax(NSpectr, 8);
    vec = [vec; nonzeros(NSpectr(BW))];
end

v_center1 = vec - c_ln;
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
xlabel('$|V_{\varepsilon}^g|^2$', 'interpreter', 'latex', 'FontSize', 20);
