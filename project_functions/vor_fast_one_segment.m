function [ESeg] = vor_fast_one_segment(E_TFR, V_, R_, op)
%VOR_FAST_OSM voronoi fast one segment mean
% op is an operation on a positiv vector

[Nfft, L] = size(E_TFR);

LR = length(R_);


%% while out of bounds
oob = 1==1;
while oob > 0
    %% random selection
    id_R = randi(LR, 1);
    VPs = R_{id_R};
    L_ = length(VPs);
    r1 = randi(L_, 1);
    r2 = mod(r1, L_) + 1;

    %% oob check
    for r=[r1, r2]
        y_min = min(V_(VPs(r), :));
        n_max = max(V_(VPs(r), 1));
        k_max = max(V_(VPs(r), 2));

        if (y_min < 1) ||...
                (n_max > L) ||...
                (k_max > Nfft)
            oob = 1;
            break;
        else
            oob = 0;
        end
    end
end

v1 = VPs(r1);
v2 = VPs(r2);
u0 = V_(v1, :);
u1 = V_(v2, :);
[bx, by] = bresenham(u0(1), u0(2), u1(1), u1(2));
b_len = length(bx);

e_vec = zeros(1, b_len);
for p=1:b_len
    e_vec(p) = E_TFR(by(p), bx(p));
end

ESeg = op(e_vec);

end

