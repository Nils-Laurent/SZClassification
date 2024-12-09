function [t, Nfft, sigs, s_sigma, s_names, VT_all] = get_signals(snr_in)
% VT_all : 0 noise, 1 signal, 2 interf;
Nfft = 1024;
L = 1024;
t = (0:(L - 1))/L;

sigs = zeros(3, L);
VT_all = zeros(3, Nfft, L);
id = 0;

% chi2inv(0.999, 2) = 13.8155;
% Y_inf = sqrt(13.8155);
% Y_sup = sqrt(13.8155);
Y_inf = 3;
Y_sup = 3;

    function [V_Type] = internal_reg(i_s1, i_s2, i_g, i_inf, i_sup)
        STFT1 = stft(i_s1, Nfft, i_g);
        STFT2 = stft(i_s2, Nfft, i_g);
        R1 = (abs(STFT1) > i_inf);
        R2 = (abs(STFT2) > i_inf);
        V_Type = R1 + 2*R2;
    end

%% signal LC
id = id + 1;
phi1 = 180*t + 250*t.^2;
phi2 = 220*t + 250*t.^2;
s1 = exp(2i*pi*phi1);
s2 = exp(2i*pi*phi2);

B = 500;
s_sigma(id) = 1/sqrt(B);

s_LC = s1 + s2;
sigs(id, :) = s_LC;
s_names(id) = "LC";

g = gauss_win(L, s_sigma(id));
nL2g = sqrt(sum(g.^2));

[~, ~, n_inf] = add_noise(s1 + s2, snr_in);
[~, ~, n_sup] = add_noise(s1 + s2, snr_in);

c_inf = Y_inf*n_inf*nL2g;
c_sup = Y_sup*n_sup*nL2g;
% gamma = n_inf*nL2g
% c_sup - c_inf
VT_all(id, :, :) = internal_reg(s1, s2, g, c_inf, c_sup);

%% signal LC cross
id = id + 1;
phi1 = 430*t + (120)/2*t.^2;
phi2 = 500*t;
s1 = exp(2i*pi*phi1);
s2 = exp(2i*pi*phi2);

B = 500;
s_sigma(id) = 1/sqrt(B);

s_LC_cross = s1 + s2;
sigs(id, :) = s_LC_cross;
s_names(id) = "crossing_LC";

g = gauss_win(L, s_sigma(id));
nL2g = sqrt(sum(g.^2));

[~, ~, n_inf] = add_noise(s1 + s2, snr_in);
[~, ~, n_sup] = add_noise(s1 + s2, snr_in);

c_inf = Y_inf*n_inf*nL2g;
c_sup = Y_sup*n_sup*nL2g;
VT_all(id, :, :) = internal_reg(s1, s2, g, c_inf, c_sup);

%% signal Oscillatory
id = id + 1;
phi1 = 475*t + 40*cos(2*pi*t);
phi2 = 525*t + 40*cos(2*pi*t);
s1 = exp(2i*pi*phi1);
s2 = exp(2i*pi*phi2);

B = 500;
s_sigma(id) = 1/sqrt(B);

s_osc = s1 + s2;
sigs(id, :) = s_osc;
s_names(id) = "osc";

g = gauss_win(L, s_sigma(id));
nL2g = sqrt(sum(g.^2));

[~, ~, n_inf] = add_noise(s1 + s2, snr_in);
[~, ~, n_sup] = add_noise(s1 + s2, snr_in);

c_inf = Y_inf*n_inf*nL2g;
c_sup = Y_sup*n_sup*nL2g;
VT_all(id, :, :) = internal_reg(s1, s2, g, c_inf, c_sup);
end

