function [id_int_wrong] = vor_find_wrong_int(id_int, id_n, R, nv, kv)

LR = length(R);
L_int = length(id_int);

id_int_wrong = zeros(1, LR);
id_d = zeros(1, LR);

N_TABLE = zeros(L_int, LR);

for z=1:LR
    verts_r = R{z};
    for j_int = 1:L_int
        z_int = id_int(j_int);
        
        if z == z_int
            continue;
        end
        
        verts_int = R{z_int};
        if intersect(verts_r, verts_int)
            N_TABLE(j_int, z) = z;
        end
    end
end

for j_int=1:L_int
    z_int = id_int(j_int);
    
%     vec_z = nonzeros(N_TABLE(j_int, :));
    vec_z_n = nonzeros(N_TABLE(j_int, id_n));
    
    if ~isempty(vec_z_n)
        id_int_wrong(z_int) = z_int;
    end
%     figure;
%     voronoi(nv, kv);
%     hold on;
%     plot(nv(z_int), kv(z_int), 'co');
%     plot(nv(vec_z), kv(vec_z), 'gx');
%     plot(nv(vec_z_n), kv(vec_z_n), 'ro');
%     hold off;
end

id_int_wrong = nonzeros(id_int_wrong);
% id_d = nonzeros(id_d);
end

