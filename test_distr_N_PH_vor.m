close all;

addpath('./TF-Toolbox/TF_Toolbox');
addpath('./TF-Toolbox/TF_Toolbox/fig/');

%% init
Nfft = 1024;
L = 1024;
t = (0:(L - 1))/L;
Fx = (0:(L-1))*L/Nfft;

%% signal def
phi1 = 490*t;
phi2 = 510*t;
s1 = exp(2i*pi*phi1);
s2 = exp(2i*pi*phi2);

B = 500;
sigma_w = 1/sqrt(B);

s = s1 + s2;

SNR_init = 5;

g = gauss_win(L, sigma_w);
g_L2Norm = sqrt(sum(g.^2));

TFRsc(stft(s_noise, Nfft, g));
% pause;

%% compute energy at cell borders
% N_rep = 1000;
N_rep = 100;
r = 1;
rr = 1;

% while r <= -1
vec = [];
while r <= N_rep
    if mod(r, 10) == 0
        fprintf("%u/%u\n", r, N_rep);
    end
    
%     nri_std = 1/sqrt(2);
%     n_var = 2*(nri_std^2);
%     
%     noise = nri_std*randn(1, L) + nri_std*1i*randn(1, L);
    s_noise = add_noise(s, SNR_init);
    
    STFT = stft(s_noise, Nfft, g);
    gamma_hat = noise_level(STFT);
    NSpectr = abs(STFT).^2/gamma_hat^2;
    
    try
        [EVM, RE] = vor_segment_operation(NSpectr, @(x)(max(x)));
    catch ME
        continue;
    end
    r = r + 1;
    rr = rr + 1;
    vec = [vec; nonzeros(EVM)];
    
%     if mod(rr, 400) == 0
%         fprintf("pause 60 sec\n");
%         pause(60);
%     end
end

% save('data_vor_max.mat', 'vec');
% load('data_vor_max.mat');
% load('data_vor_max_1000.mat');

yMin = min(vec);
yMax = max(vec);
dx = (yMax - yMin)/500;
x = yMin:dx:yMax;

%% fit EV
% [p_hat, p_ic] = evfit(vec);
% p_mu = p_hat(1);
% p_sig = p_hat(2);
% 
% g_pdf = evpdf(vec, p_mu, p_sig);

%% fit GEV
% [p_hat, p_ic] = gevfit(vec);
% 
% p_hat = [shape, scale, location]
% %  When p_K > 0, the GEV distribution is the type II, or Frechet, extreme value distribution.
% p_K = p_hat(1);
% p_sig = p_hat(2);
% p_mu = p_hat(3);
% 
% g_pdf = gevpdf(x, p_K, p_sig, p_mu);

%% Weibull
% DObj = fitdist(vec, 'Weibull');
% g_pdf = DObj.pdf(x);

%% Gamma
% DObj = fitdist(vec, 'gamma');
p_shape = 2.2792;
p_scale = 1.7074;
DObj = makedist('gamma', 'a', p_shape, 'b', p_scale);
g_pdf = DObj.pdf(x);

DObj.ParameterDescription
DObj.ParameterNames
DObj.ParameterValues

%% figures
figure;
histogram(vec, 'normalization', 'pdf', 'BinWidth', 0.25);
xlabel('$|V_{\varepsilon}^g|^2$', 'interpreter', 'latex', 'FontSize', 20);
hold on;
plot(x, g_pdf, 'r-');
hold off;