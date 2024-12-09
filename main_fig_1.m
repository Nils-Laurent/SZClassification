%% Distribution of zeros and local max
close all;

addpath('./TF-Toolbox/TF_Toolbox');
addpath('./TF-Toolbox/TF_Toolbox/fig/');
addpath('./project_functions');

%% init
Nfft = 1024;
L = 1024;
Tx = (0:(L - 1))/L;
Fx = (0:(L-1))*L/Nfft;

B = 500;
sigma_w = 1/sqrt(B);
g = gauss_win(L, sigma_w);
g_L2Norm = sqrt(sum(g.^2));

n_std = 1;
noise = n_std*randn(1, L) + n_std*1i*randn(1, L);

STFT = stft(noise, Nfft, g);
NSpectr = abs(STFT).^2/(n_std*g_L2Norm)^2;

[nv, kv] = zeros_spec(NSpectr);
[nv2, kv2] = local_max_spec(NSpectr);

TFRsc(Tx, Fx, NSpectr)
legend('location', 'northeast');
hold on;
plot(Tx(nv), Fx(kv), 'bo', 'DisplayName', 'Zeros');
plot(Tx(nv2), Fx(kv2), 'rx', 'DisplayName', 'Local maxima');
hold off;
xlim([0, 0.5]);
ylim([0, 512]);
set(gca,'XTick',[], 'YTick', []);
write_figfiles('FIG1_local_max');

%% compute energy at cell borders
N_rep = 300;

vec = [];
for r = 1:N_rep
    if mod(r, 10) == 0
        fprintf("%u/%u\n", r, N_rep);
    end
    
    noise = randn(1, L) + 1i*randn(1, L);
    STFT = stft(noise, Nfft, g);
    NSpectr = abs(STFT).^2/(g_L2Norm^2);
    
    BW = imregionalmax(NSpectr, 8);
    vec = [vec; nonzeros(NSpectr(BW))];
end



yMin = min(vec);
yMax = max(vec);
dx = (yMax - yMin)/500;
x = yMin:dx:yMax;

%% fit EV (gumbel)
GumbelObj = fitdist(vec, 'ev');
Gumbel_pdf = GumbelObj.pdf(x);

fprintf("Gumbel\n");
for id=1:length(GumbelObj.ParameterDescription)
    fprintf("%s(%s) = %f\n", GumbelObj.ParameterNames{id},...
        GumbelObj.ParameterDescription{id},...
        GumbelObj.ParameterValues(id));
end

%% Gamma
GammaObj = fitdist(vec, 'gamma');
Gamma_pdf = GammaObj.pdf(x);

fprintf("Gamma\n");
for id=1:length(GammaObj.ParameterDescription)
    fprintf("%s(%s) = %f\n", GammaObj.ParameterNames{id},...
        GammaObj.ParameterDescription{id},...
        GammaObj.ParameterValues(id));
end

%% figures
fig_form;
histogram(vec, 'normalization', 'pdf', 'DisplayName', 'Local max energy');
xlabel('$|V_{\varepsilon}^g|^2$');
hold on;
plot(x, Gumbel_pdf, 'b--', 'DisplayName', 'Gumbel PDF');
plot(x, Gamma_pdf, 'r-', 'DisplayName', 'Gamma PDF');
hold off;
legend('location', 'northeast');
write_figfiles("FIG1_LM_gamma_300R");