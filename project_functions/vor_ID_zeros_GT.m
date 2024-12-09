function [id_zero_GT] = vor_ID_zeros_GT(GT_TFR, nv, kv)

    id_zero_GT = zeros(1, length(nv));

    GT_sig1 = (GT_TFR == 1) | (GT_TFR == 3);
    [~, Re_R1, ~, R_struct, id_border] = ...
        vor_segment_operation(GT_sig1, @(x)(max(x)), nv, kv);
    
    GT_sig2 = GT_TFR >= 2;
    [~, Re_R2, ~, ~, ~] = vor_segment_operation(GT_sig2, @(x)(max(x)), nv, kv);


    for j=1:length(Re_R1)
        if id_border(j) > 0
            id_zero_GT(j) = -1;
            continue;
        end
        
        Ej_R1 = Re_R1{j};
        Ej_R2 = Re_R2{j};

        ctr_R1 = 0;
        ctr_R2 = 0;
        ctr_any = 0;
        for r = 1:length(Ej_R1)
            ctr_R1 = ctr_R1 + (Ej_R1(r) > 0);
            ctr_R2 = ctr_R2 + (Ej_R2(r) > 0);
            ctr_any = ctr_any + ((Ej_R1(r) > 0) || (Ej_R2(r) > 0));
        end

        if ctr_R1 > 0 && ctr_R2 > 0
%         if (ctr_any == length(Ej_R1)) && ...
%                 (ctr_R1 > 0 && ctr_R2 > 0)
            id_zero_GT(j) = 2;
        elseif ctr_R1 > 2 || ctr_R2 > 2
            id_zero_GT(j) = 1;
        else
            id_zero_GT(j) = 0;
        end
    end
    
    id_int = find(id_zero_GT == 2);
    id_sig = find(id_zero_GT == 1);
    id_n = find(id_zero_GT == 0);
    
    [id_int_wrong] = vor_find_wrong_int(id_int, id_n, R_struct, nv, kv);
    
    id_zero_GT(id_int_wrong) = 1;
    
%     TFRsc(GT_TFR);
%     hold on;
%     voronoi(nv, kv);
%     plot(nv(id_sig), kv(id_sig), 'rx');
%     plot(nv(id_int), kv(id_int), 'co');
%     plot(nv(id_int_wrong), kv(id_int_wrong), 'go');
%     hold off;
end

