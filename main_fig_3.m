%% Test of the bias

close all;

addpath('./TF-Toolbox/TF_Toolbox');
addpath('./TF-Toolbox/TF_Toolbox/fig/');
addpath('./project_functions');

%% test
[t, Nfft, sigs, s_sigma, s_names, ~] = get_signals(0);

L = length(t);
freq = (0:(L-1))*L/Nfft;

for id = 1:length(s_sigma)
    s_names(id)
    [SNRs, p_gh, p_g_trunc] = main_fig_3_test(sigs(id, :), s_sigma(id), Nfft);
    
    %% figure accuracy
    fig_form;
    hold on;
    plot(SNRs, p_gh, '--', 'DisplayName', 'Error of $\widehat{\gamma}$');
    plot(SNRs, p_g_trunc, '-', 'DisplayName', 'Error of $\tilde{\gamma}$');
    hold off;
    xlabel("input SNR");
    ylabel("estimation error");
    legend('location', 'northwest');
    write_figfiles(['FIG3_est_error_50R_', s_names(id)]);
    close all;
end

%% figures



% yMin = min(vec);
% yMax = max(vec);
% dx = (yMax - yMin)/500;
% x = yMin:dx:yMax;
% 
% g_pdf = DObj.pdf(x);
% 
% figure;
% histogram(vec, 'normalization', 'pdf', 'BinWidth', 0.25);
% xlabel('$|V_{\varepsilon}^g|^2$', 'interpreter', 'latex', 'FontSize', 20);
% hold on;
% plot(x, g_pdf, 'r-');
% hold off;

% figure;
% histogram(vec2, 'normalization', 'pdf', 'BinWidth', 0.25);
% xlabel('$|V_{\varepsilon}^g|^2$', 'interpreter', 'latex', 'FontSize', 20);
% hold on;
% plot(x, g_pdf, 'g-');
% hold off;