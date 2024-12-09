%% Distribution of energy on voronoi edges in PURE NOISE

close all;

addpath('./TF-Toolbox/TF_Toolbox');
addpath('./TF-Toolbox/TF_Toolbox/fig/');
addpath('./project_functions');

%% init
Nfft = 1024;
L = 1024;
t = (0:(L - 1))/L;
Fx = (0:(L-1))*L/Nfft;

B = 500;
sigma_w = 1/sqrt(B);
g = gauss_win(L, sigma_w);
g_L2Norm = sqrt(sum(g.^2));

%% compute energy at cell borders
N_rep = 300;
r = 1;

vec = [];
while r <= N_rep
    if mod(r, 10) == 0
        fprintf("%u/%u\n", r, N_rep);
    end
    
    n_std = 1;
    
    noise = n_std*randn(1, L) + n_std*1i*randn(1, L);
    
    STFT = stft(noise, Nfft, g);
    NSpectr = abs(STFT).^2/(n_std*g_L2Norm)^2;
    
    try
        [EVM, RE] = vor_segment_operation(NSpectr, @(x)(max(x)));
    catch ME
        continue;
    end
    r = r + 1;
    vec = [vec; nonzeros(EVM)];
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
histogram(vec, 'normalization', 'pdf', 'DisplayName', 'Edges energy');
xlabel('$|V_{\varepsilon}^g|^2$');
hold on;
plot(x, Gumbel_pdf, 'b--', 'DisplayName', 'Gumbel PDF');
plot(x, Gamma_pdf, 'r-', 'DisplayName', 'Gamma PDF');
hold off;
legend('location', 'northeast');
write_figfiles("FIG2_distribution_noise_300R");