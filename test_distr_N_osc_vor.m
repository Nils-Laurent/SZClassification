close all;

addpath('./TF-Toolbox/TF_Toolbox');
addpath('./TF-Toolbox/TF_Toolbox/fig/');

%% init
Nfft = 1024;
L = 1024;
t = (0:(L - 1))/L;
Fx = (0:(L-1))*L/Nfft;

%% signal def
phi1 = 485*t + 40*cos(2*pi*t);
phi2 = 515*t + 40*cos(2*pi*t);
s1 = exp(2i*pi*phi1);
s2 = exp(2i*pi*phi2);

B = 500;
sigma_w = 1/sqrt(B);

s = s1 + s2;

SNR_init = 5;

g = gauss_win(L, sigma_w);
g_L2Norm = sqrt(sum(g.^2));

% return;
%% compute energy at cell borders
N_rep = 1;
% N_rep = 1;
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
    [s_noise, ~, nscale] = add_noise(s, SNR_init);
%     s_noise = s_noise - s.';
    
    STFT = stft(s_noise, Nfft, g);
    gamma_hat = noise_level(STFT);
    gamma_GT = nscale*g_L2Norm;
    c_err = (gamma_GT/gamma_hat)^2;
    NSpectr = c_err*abs(STFT).^2/(gamma_GT^2);
    
    mask = ones(size(STFT));
    
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

% vec = c_err*vec;

%% Gamma
% DObj = fitdist(vec, 'gamma');
p_shape = 2.2792;% k
p_scale = 1.7074;% theta

DObj = makedist('gamma', 'a', p_shape, 'b', p_scale);

yMin = min(vec);
yMax = max(vec);
dx = (yMax - yMin)/500;
x = yMin:dx:yMax;

g_pdf = DObj.pdf(x);

%% Maximum likelihood scale (theta)
vec_tr = vec(vec < 3*gamma_hat);
T0 = sum(vec_tr)/(p_shape*length(vec_tr));
T_err = T0/p_scale;

gh2 = gamma_hat/sqrt(1/T_err);
fprintf("gamma_hat = %f, gamma_GT = %f, gamma_hat2 = %f\n",...
    gamma_hat, gamma_GT, gh2);


vec2 = vec*(gamma_hat/gh2)^2;
%     c_err = (gamma_GT/gamma_hat)^2;

%% figures

figure;
histogram(vec, 'normalization', 'pdf', 'BinWidth', 0.25);
xlabel('$|V_{\varepsilon}^g|^2$', 'interpreter', 'latex', 'FontSize', 20);
hold on;
plot(x, g_pdf, 'r-');
hold off;

figure;
histogram(vec2, 'normalization', 'pdf', 'BinWidth', 0.25);
xlabel('$|V_{\varepsilon}^g|^2$', 'interpreter', 'latex', 'FontSize', 20);
hold on;
plot(x, g_pdf, 'g-');
hold off;