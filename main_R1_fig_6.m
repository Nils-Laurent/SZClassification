%% SIGNAL SEPARATION

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
p_vec = 0:20:160;
% p_vec = 0:80:160;

N_rep = 30;
% N_rep = 2;

id_vec = 1:2;
snr_vec = [0, 5, 15];
% snr_vec = [5, 15];

N_id = length(id_vec);
N_snr = length(snr_vec);
N_p = length(p_vec);
% res = zeros(N_id, N_snr, N_p, N_rep);

for id = id_vec
    break
    fprintf("id = %u\n", id);
    for id_p = 1:length(p_vec)
        p = p_vec(id_p);
        fprintf("%u/%u ", id_p, length(p_vec));
        for id_snr=1:length(snr_vec)
            c_snr = snr_vec(id_snr);
            int_acc = 0;

            for nr = 1:N_rep
                [t, Nfft, sigs, s_sigma, s_names, VT_all] = get_signals_var(c_snr, p);
            
                [s_noise, ~, ~] = add_noise(sigs(id, :), c_snr);
                [g, Lh] = gauss_win(L, s_sigma(id));
                STFT = stft(s_noise, Nfft, g);
                STFT_ref = stft(sigs(id, :), Nfft, g);

%                 TFRsc(STFT);
%                 break
                
                [id_zeros,nv,kv]...
                    = vor_ID_zeros(s_noise, s_sigma(id), Nfft);
        
                res(id, id_snr, id_p, nr) = sum(kv(id_zeros == 2)) / sum(kv(id_zeros == 1));
%                 res(id, id_snr, id_p, nr) = (sum(kv(id_zeros == 2)) > 0);
            end

%             res(id, id_snr, id_p) = int_acc/n_rep;
        end
    end
    fprintf("\nid = %u [end]\n", id);
end

% save("data_R1_FIG6_v3.mat", "res");

fig_form;
p_form1 = ["-", "--", "-."];
p_form2 = ["-o", "--o", "-.o"];
p_name = ["Parallel linear chirps", "Crossing linear chirps",...
    "Oscillating phases"];

hold on;
for id_snr=1:length(snr_vec)
    res_snr = squeeze(res(:, id_snr, :, :));
    res_mean = mean(res_snr, 3);
    c_snr = snr_vec(id_snr);
    c_snr
%     res_mean

    d_plot = res_mean;

    id = 1;
    plot(p_vec, d_plot(id, :), p_form1(id_snr), 'DisplayName', p_name(id) + " " + string(c_snr) + "dB");

    id = 2;
    plot(p_vec, d_plot(id, :), p_form2(id_snr), 'DisplayName', p_name(id) + " " + string(c_snr) + "dB");
end
hold off;
xlim([p_vec(1), p_vec(end)]);
%     ylim([0, 1]);
xlabel("$\theta$");
ylabel("\#SS/\#SN");
legend('location', "northeast");
set(gcf,'Position',[300 300 900 600])
write_figfiles('FIG6_separ_' + string(N_rep) + 'R');

return;

%% old figure
for id_snr=1:length(snr_vec)
    res_snr = squeeze(res(:, id_snr, :, :));
    res_mean = mean(res_snr, 3);
    c_snr = snr_vec(id_snr);
    c_snr
%     res_mean

    d_plot = 2*res_mean;

    p_form = ["-", "--", "-."];
    p_form_r = ["-o", "--o", "-.o"];
    p_name = ["Parallel linear chirps", "Crossing linear chirps",...
        "Oscillating phases"];
    fig_form;
    hold on;
    for id=1:length(s_sigma)
        plot(p_vec, d_plot(id, :), p_form(id), 'DisplayName', p_name(id));
    end
    hold off;
    xlim([p_vec(1), p_vec(end)]);
%     ylim([0, 1]);
    xlabel("$\theta$");
    ylabel("ratio");
    legend('location', "northeast");
    set(gcf,'Position',[300 300 900 450])
    write_figfiles('FIG6_separ_' + string(N_rep) + 'R_' + string(c_snr) + 'SNR');
end
