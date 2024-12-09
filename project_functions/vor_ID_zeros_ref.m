function [id_zero,nv,kv] = vor_ID_zeros_ref(s_noise, sigma_w, Nfft, n_scale)

    L = length(s_noise);

    [g, Lh] = gauss_win(L, sigma_w);
    g_L2Norm = sqrt(sum(g.^2));
    
    gamma_GT = n_scale*g_L2Norm;
    gamma_est = gamma_GT;

    STFT = stft(s_noise, Nfft, g);
    NSpectr = abs(STFT).^2/(gamma_est^2);
    
    [nv, kv] = zeros_spec(NSpectr);


    %% classify cells
    [~, Re_Max, ~, R_Max, id_border] = ...
        vor_segment_operation(NSpectr, @(x)(max(x)), nv, kv);

    % GAMMA THRESHOLD : see test_distr_vor_segment.m
    [pmax_shape, pmax_scale] = vor_energy_d_param();

    %% suite

    p0 = 0.999;
    gx_max_zero = fzero(@(x)(cdf('gamma', x, pmax_shape, pmax_scale)) - p0, 34);

    % gx_zero
    TH_emax = gx_max_zero;%*(gamma_est)^2
    
%     TFRsc(abs(STFT) > 3*gamma_est);
%     hold on;
%     voronoi(nv, kv);
%     hold off;
%     
%     TFRsc(NSpectr > TH_emax);
%     hold on;
%     voronoi(nv, kv);
%     hold off;

    id_zero = zeros(1, length(nv));

    for j=1:length(R_Max)
        if id_border(j) > 0
            id_zero(j) = -1;
            continue;
        end
        
        Ej_max = Re_Max{j};

        %% counter of high energy segments
        ctr_MAX = 0;
        for r = 1:length(Ej_max)
            ctr_MAX = ctr_MAX + (Ej_max(r) > TH_emax);
        end

        if ctr_MAX == length(Ej_max)
            id_zero(j) = 2;
        elseif ctr_MAX > 2
            id_zero(j) = 1;
        end
    end
    
    id_int = find(id_zero == 2);
    id_n = find(id_zero == 0);
    
    [id_int_wrong] = vor_find_wrong_int(id_int, id_n, R_Max, nv, kv);
    
    id_zero(id_int_wrong) = 1;
end

