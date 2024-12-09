function [M] = zero_mask(TFR)

SA = abs(TFR);
M = imregionalmin(SA, 8);

end

