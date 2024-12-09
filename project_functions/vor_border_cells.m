function [id_border, id_inner] = vor_border_cells(V_coord, R, L, Nfft)

LR = length(R);

id_border = zeros(1, LR);

for j=1:LR
    V_id_vec = R{j};
    
    for r = 1:(length(V_id_vec))
        vr = V_coord(V_id_vec(r), :);
        
        if vr(1) < 1 ...
            || vr(1) > L ...
            || vr(2) < 1 ...
            || vr(2) > Nfft
        
            id_border(j) = 1;
            break;
        end
    end
end

id_border = id_border > 0;
id_inner = id_border == 0;

end

