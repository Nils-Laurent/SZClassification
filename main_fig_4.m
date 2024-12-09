%% Test of the bias

close all;

addpath('./TF-Toolbox/TF_Toolbox');
addpath('./TF-Toolbox/TF_Toolbox/fig/');
addpath('./project_functions');

%% test
c_snr = 8;
[t, Nfft, sigs, s_sigma, s_names, VT_all] = get_signals(c_snr);
L = length(sigs(1, :));
freq = (0:Nfft)*L/Nfft;

for id=1:length(s_sigma)
%     break;
    [s_noise, ~, ~] = add_noise(sigs(id, :), c_snr);
    [g, Lh] = gauss_win(L, s_sigma(id));
    STFT = stft(s_noise, Nfft, g);
    STFT_ref = stft(sigs(id, :), Nfft, g);
    
    [id_zeros,nv,kv]...
        = vor_ID_zeros(s_noise, s_sigma(id), Nfft);

    id_sig = id_zeros == 1;
    id_int = id_zeros == 2;
    
    [id_zero_GT] = vor_ID_zeros_GT(squeeze(VT_all(id, :, :)), nv, kv);
    
    
    L = length(s_noise);
    [g, Lh] = gauss_win(L, s_sigma(id));
    
    nv_sigGT = nv(id_zero_GT == 1);
    kv_sigGT = kv(id_zero_GT == 1);
    nv_intGT = nv(id_zero_GT == 2);
    kv_intGT = kv(id_zero_GT == 2);
    
    Tx=1:L;
    Fx=1:Nfft;
    
    TFRsc(Tx, Fx, abs(STFT_ref));
%     xlim([Lh, L - Lh]);
    set(gca,'XTick',[], 'YTick', []);
    write_figfiles(['FIG4_TFR_', s_names(id)]);

    TFRsc(Tx, Fx, abs(STFT));
    hold on;
    voronoi(nv, kv);
    plot(nv(id_sig), kv(id_sig), 'rx');
    plot(nv(id_int), kv(id_int), 'co');
    hold off;
    xlim([Lh, L - Lh]);
    set(gca,'XTick',[], 'YTick', []);
    write_figfiles(['FIG4_illu_', s_names(id)]);
    
%     TFRsc(squeeze(VT_all(id, :, :)));
%     hold on;
%     voronoi(nv, kv);
%     plot(nv_sigGT, kv_sigGT, 'rx');
%     plot(nv_intGT, kv_intGT, 'co');
%     hold off;
% %     close all;
end

% return;
ev_all = [];
ev_all_r = [];
dNN_all = [];
dSN_all = [];
dSS_all = [];
for id=1:length(s_sigma)
    [SNRs, err_vec, err_vec_r, dNN, dSN, dSS] =...
        main_fig_4_test(sigs(id, :), s_sigma(id), Nfft, id);
    ev_all(id, :) = err_vec;
    ev_all_r(id, :) = err_vec_r;
    dNN_all(id, :) = dNN;
    dSN_all(id, :) = dSN;
    dSS_all(id, :) = dSS;
end
    
p_form = ["-", "--", "-."];
p_form_r = ["-o", "--o", "-.o"];
p_name = ["Parallel linear chirps", "Crossing linear chirps",...
    "Oscillating phases"];
fig_form;
hold on;
for id=1:length(s_sigma)
    plot(SNRs, 100*ev_all(id, :), p_form(id), 'DisplayName', p_name(id));
end
hold off;
xlabel("input SNR");
ylabel("classification error (\%)");
legend('location', "northwest");
write_figfiles('FIG4_valid_30R');

return;

%% old figures
fig_form;
hold on;
id = 1;
plot(SNRs, dNN_all(id, :), '-', 'DisplayName', 'dNN');
plot(SNRs, dSN_all(id, :), '--', 'DisplayName', 'dSN');
plot(SNRs, dSS_all(id, :), '-.', 'DisplayName', 'dSS');
hold off;
xlabel("input SNR");
ylabel("estimation error");
legend('location', "northwest");
write_figfiles('FIG4_valid_30R_id1');

fig_form;
hold on;
id = 2;
plot(SNRs, dNN_all(id, :), '-', 'DisplayName', 'dNN');
plot(SNRs, dSN_all(id, :), '--', 'DisplayName', 'dSN');
plot(SNRs, dSS_all(id, :), '-.', 'DisplayName', 'dSS');
hold off;
xlabel("input SNR");
ylabel("estimation error");
legend('location', "northwest");
write_figfiles('FIG4_valid_30R_id2');

fig_form;
hold on;
id = 3;
plot(SNRs, dNN_all(id, :), '-', 'DisplayName', 'dNN');
plot(SNRs, dSN_all(id, :), '--', 'DisplayName', 'dSN');
plot(SNRs, dSS_all(id, :), '-.', 'DisplayName', 'dSS');
hold off;
xlabel("input SNR");
ylabel("estimation error");
legend('location', "northwest");
write_figfiles('FIG4_valid_30R_id3');
