close all;

addpath('./TF-Toolbox/TF_Toolbox');
addpath('./TF-Toolbox/TF_Toolbox/fig/');

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
s_noise = add_noise(s, SNR_init);

%% time frequency

g = gauss_win(L, sigma_w);
g_L2Norm = sqrt(sum(g.^2));

STFT = stft(s_noise, Nfft, g);
gamma_hat = noise_level(STFT);


% cr_std = 1/sqrt(2);
% ci_std = cr_std;
% n_var = cr_std.^2 + ci_std.^2;
% GT = sqrt(n_var)/sqrt(2)
% 
% NR = 10;
% g_rep = 0;
% for r = 1:NR
%     if mod(r, 10) == 0
%         fprintf("%u/%u, ", r, NR);
%     end
%     noise = cr_std*randn(1, L) + ci_std*1i*randn(1, L);
%     STFT = stft(noise, Nfft, g);
%     gamma_hat = noise_level(STFT);
%     g_rep = g_rep + gamma_hat;
% end
% fprintf("\n");
% 
% g_rep = g_rep/NR;
% g_rep = g_rep/g_L2Norm;
% 
% Est1 = g_rep
% 
% return;
%% classify cells
spectro = abs(STFT).^2;

[nv, kv] = zeros_spec(STFT);

DT = delaunayTriangulation(nv, kv);
[V_, R_] = voronoiDiagram(DT);
LR = length(R_);

[EVM, R_e] = vor_segment_operation(spectro, @(x)(max(x)));
vec = [vec; nonzeros(EVM)];

% BLIND THRESHOLD
% TH_e = (4*noise_level(STFT))^2;

% GAMMA THRESHOLD : see test_distr_vor_segment.m
p_shape = 2.2792;
p_scale = 1.7074;
DObj = makedist('gamma', 'a', p_shape, 'b', p_scale);

%% Maximum likelihood scale (theta)
vec_tr = vec(vec < 3*gamma_hat);
T0 = sum(vec_tr)/(p_shape*length(vec_tr));
T_err = T0/p_scale;

gh2 = gamma_hat/sqrt(1/T_err);
fprintf("gamma_hat = %f, gamma_GT = %f, gamma_hat2 = %f\n",...
    gamma_hat, gamma_GT, gh2);


vec2 = vec*(gamma_hat/gh2)^2;

%% figures
figure;
histogram(vec, 'normalization', 'pdf');
xlabel('$|V_{\varepsilon}^g|^2$', 'interpreter', 'latex', 'FontSize', 20);
hold on;
plot(x, g_pdf, 'r-');
hold off;

%% suite

p0 = 0.99;
gx_zero = fzero(@(x)(cdf('gamma', x, p_shape, p_scale)) - p0, 34);

% gx_zero
TH_e = gx_zero*(gamma_hat)^2

nv_sig = [];
kv_sig = [];

nv_int = [];
kv_int = [];

for j=1:LR
    Ej = R_e{j};
    
    %% counter of high energy segments
    ctr = 0;
    for r = 1:length(Ej)
        b_r = Ej(r) > TH_e;
        
        if b_r > 0
            ctr = ctr + 1;
        else
            ctr = 0;
        end
    end
    
    if ctr == length(Ej)
        nv_int = [nv_int, nv(j)];
        kv_int = [kv_int, kv(j)];
    elseif ctr > 0
        nv_sig = [nv_sig, nv(j)];
        kv_sig = [kv_sig, kv(j)];
    end
end

TFRsc(spectro);
hold on;
voronoi(nv, kv);
plot(nv_sig, kv_sig, 'rx');
plot(nv_int, kv_int, 'co');
hold off;

% NClass = 3;
% IDX = kmeans(EReg.', NClass, 'Replicates', 10);
