function [SNRs, p_ghat, p_g_trunc] = main_fig_3_test(s_ref, sigma_w, Nfft)
L = length(s_ref);
g = gauss_win(L, sigma_w);
g_L2Norm = sqrt(sum(g.^2));

%% compute energy at cell borders
SNRs = 0:2:20;
% SNRs = 5;
N_snr = length(SNRs);
p_ghat = zeros(1, N_snr);
p_g_trunc = zeros(1, N_snr);

N_rep = 50;

for n_snr = 1:N_snr
    c_snr = SNRs(n_snr);
    fprintf("snr %u/%u\n", n_snr, N_snr);
    
    r = 1;
    while r <= N_rep
        if mod(r, 10) == 0
            fprintf("rep %u/%u\n", r, N_rep);
        end

        [s_noise, ~, nscale] = add_noise(s_ref, c_snr);

        STFT = stft(s_noise, Nfft, g);
        gamma_hat = noise_level(STFT);
        gamma_GT = nscale*g_L2Norm;
        
        STFT_tr = STFT(abs(STFT) < 3*gamma_hat);
        g_trunc = noise_level(STFT_tr);
        
        p_ghat(n_snr) = p_ghat(n_snr) + abs(gamma_GT - gamma_hat)/gamma_GT;
        p_g_trunc(n_snr) = p_g_trunc(n_snr) + abs(gamma_GT - g_trunc)/gamma_GT;
        
        r = r + 1;
    end
    
    p_ghat(n_snr) = p_ghat(n_snr)/N_rep;
    p_g_trunc(n_snr) = p_g_trunc(n_snr)/N_rep;
end

end


