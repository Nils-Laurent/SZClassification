function [nv, kv] = zeros_spec(TFR)

M = zero_mask(TFR);
[N, L] = size(TFR);

nv = zeros(size(TFR));
kv = zeros(size(TFR));
for n=2:(L - 1)
    for k=2:(N - 1)
        if M(k, n) > 0
            nv(k, n) = n;
            kv(k, n) = k;
        end
    end
end

nv = nonzeros(nv(:));
kv = nonzeros(kv(:));

end

