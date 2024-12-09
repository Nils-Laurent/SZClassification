function [SNRs, vec_err, vec_err_ref, dNN, dSN, dSS] = main_fig_4_test(s_ref, sigma_w, Nfft, id)
% L = length(s_ref);
% g = gauss_win(L, sigma_w);

%% compute energy at cell borders
SNRs = 0:2:20;
% SNRs = 0;
N_snr = length(SNRs);

N_rep = 30;
% N_rep = 1;

vec_err = zeros(1, N_snr);
vec_err_ref = vec_err;
dNN = vec_err;
dSN = vec_err;
dSS = vec_err;

for n_snr = 1:N_snr
    c_snr = SNRs(n_snr);
    [~, ~, ~, ~, ~, VT_all] = get_signals(c_snr);
    VT = squeeze(VT_all(id, :, :));
    
    fprintf("snr %u/%u\n", n_snr, N_snr);
    
    r = 1;
    while r <= N_rep
        if mod(r, 10) == 0
            fprintf("rep %u/%u\n", r, N_rep);
        end

        [s_noise, ~, n_scale] = add_noise(s_ref, c_snr);
        
%         [id_zeros, nv, kv]...
%             = vor_ID_zeros(s_noise, sigma_w, Nfft);
        try
            [id_zeros, nv, kv]...
                = vor_ID_zeros(s_noise, sigma_w, Nfft);
        catch RE
            continue
        end
        
        %% disable exp. with true gamma
%         [id_zeros_r, ~]...
%             = vor_ID_zeros_ref(s_noise, sigma_w, Nfft, n_scale);
        

        %% remove zeros on the border
        [id_zero_GT] = vor_ID_zeros_GT(VT, nv, kv);
        
        %% ignore zeros on the border
        L = length(s_noise);
        [g, Lh] = gauss_win(L, sigma_w);

%         nv_sig = nv(id_zeros == 1);
%         kv_sig = kv(id_zeros == 1);
%         nv_int = nv(id_zeros == 2);
%         kv_int = kv(id_zeros == 2);
%         STFT = stft(s_noise, Nfft, g);
%         TFRsc(STFT);
%         hold on;
%         voronoi(nv, kv);
%         plot(nv_sig, kv_sig, 'rx');
%         plot(nv_int, kv_int, 'co');
%         hold off;
%         
%         nv_sigGT = nv(id_zero_GT == 1);
%         kv_sigGT = kv(id_zero_GT == 1);
%         nv_intGT = nv(id_zero_GT == 2);
%         kv_intGT = kv(id_zero_GT == 2);
%         TFRsc(VT);
%         hold on;
%         voronoi(nv, kv);
%         plot(nv_sigGT, kv_sigGT, 'rx');
%         plot(nv_intGT, kv_intGT, 'co');
%         hold off;
        
        id_relz = id_zeros >= 0;
        L_relz = sum(id_relz);

        vec_err(n_snr) = vec_err(n_snr) +...
            sum(id_zeros(id_relz) == id_zero_GT(id_relz))/L_relz;
        
%         vec_err_ref(n_snr) = vec_err_ref(n_snr) +...
%             sum(id_zeros_r(L_relz) == id_zero_GT(L_relz))...
%             /L_relz;
        
        i_dNN = abs(sum(id_zeros == 0) - sum(id_zero_GT == 0));
        i_dSN = abs(sum(id_zeros == 1) - sum(id_zero_GT == 1));
        i_dSS = abs(sum(id_zeros == 2) - sum(id_zero_GT == 2));
        
        dNN(n_snr) = dNN(n_snr) + i_dNN/L_relz;
        dSN(n_snr) = dSN(n_snr) + i_dSN/L_relz;
        dSS(n_snr) = dSS(n_snr) + i_dSS/L_relz;
        
        r = r + 1;
    end
end
vec_err = 1 - vec_err/N_rep;
% vec_err_ref = 1 - vec_err_ref/N_rep;
dNN = dNN/N_rep;
dSN = dSN/N_rep;
dSS = dSS/N_rep;

end


