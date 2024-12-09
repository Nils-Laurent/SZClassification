% statistique max

close all;

l_spec = 1/2;
N_rep = 4096;
N_gumbel = 4096;

v_exp_max = max(exprnd(1/l_spec, N_gumbel, N_rep));
c_ln = log(N_gumbel)/l_spec;
v_center = v_exp_max - c_ln;

%% gumbel parameters
g_sig = 1/l_spec;
g_mu = 0;

g_max = max(v_center);
g_min = min(v_center);
x = g_min:g_max/500:g_max;

%% gumbel pdf, formulation de P. Flandrin
z = -(x - g_mu)/g_sig;
g_pdf = 1/g_sig*exp(z - exp(z));

figure;
hold on;
histogram(v_center, 'Normalization', 'pdf');
plot(x, g_pdf, 'r-');
hold off;