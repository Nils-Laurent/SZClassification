function [EVM, R_energy, V_coord, R_struct, id_border] = ...
    vor_segment_operation(E_TFR, op, nv_i, kv_i)
% op is an operation on a positiv vector

[Nfft, L] = size(E_TFR);

if nargin > 2
    nv = nv_i;
    kv = kv_i;
else
    [nv, kv] = zeros_spec(E_TFR);
end

DT = delaunayTriangulation(nv(:), kv(:));
[V_coord, R_struct] = voronoiDiagram(DT);
R_energy = R_struct;

%% identify border cells
[id_border, ~] = vor_border_cells(V_coord, R_struct, L, Nfft);

LR = length(R_struct);
LV = length(V_coord);


%% compute energy at cell borders

% EVM : vertex energy matrix
%       triangular (symmetry) [(x1, y1), (x2, y2)] = [(x2, y2), (x1, y1)]
%       ensures unicity
EVM = zeros(LV);

for j=1:LR
    
    if id_border(j) > 0
        continue;
    end
    
    edge_vert_ids = R_struct{j};
    L_Rverts = length(edge_vert_ids);
    
    for r = 1:L_Rverts
        v1 = edge_vert_ids(r);
        v2 = edge_vert_ids(mod(r, L_Rverts) + 1);
        coord_v1 = V_coord(v1, :);
        coord_v2 = V_coord(v2, :);
        [bx, by] = bresenham(coord_v1(1), coord_v1(2),...
            coord_v2(1), coord_v2(2));
        b_len = length(bx);

        vec_e = zeros(1, b_len);
        for p=1:b_len
            vec_e(p) = E_TFR(by(p), bx(p));
        end
        v_e = op(vec_e);
        
        EVM(v1, v2) = v_e;
        EVM(v2, v1) = v_e;
        R_energy{j}(r) = v_e;
    end
end

for j=1:LV
    for k=j+1:LV
        EVM(k, j) = 0;
    end
end

% vor_mean = nonzeros(EVM(:));
end

