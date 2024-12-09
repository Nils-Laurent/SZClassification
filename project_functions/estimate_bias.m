function [gamma_new, T_err] = estimate_bias(vec, gamma_hat)

[p_shape, p_scale] = vor_energy_d_param();

% vec_tr = vec(vec < 3*gamma_hat);
vec_tr = vec(vec < 9);
T0 = sum(vec_tr)/(p_shape*length(vec_tr));
T_err = T0/p_scale;

% new_gamma = gamma_hat/sqrt(1/T_err);
gamma_new = gamma_hat*sqrt(T_err);

        
%% debug
% vec2 = vec/T_err;
% [p_shape, p_scale] = vor_energy_d_param();
% DObj = makedist('gamma', 'a', p_shape, 'b', p_scale);
% 
% yMin = min(vec2);
% yMax = max(vec2);
% dx = (yMax - yMin)/10000;
% x = yMin:dx:yMax;
% 
% figure;
% %         histogram(vec/T_err, 'normalization', 'pdf', 'BinWidth', 0.25);
% histogram(vec2, 'normalization', 'pdf', 'BinWidth', 0.25);
% xlabel('$|V_{\varepsilon}^g|^2$', 'interpreter', 'latex', 'FontSize', 20);
% hold on;
% plot(x, DObj.pdf(x), 'r-');
% plot([3*gamma_hat, 3*gamma_hat], [0, 0.2], 'c');
% hold off;
% xlim([0, 15]);
% pause;

end

