function [M] = LM_mask(TFR)

SA = abs(TFR);
M = imregionalmax(SA, 8);

end

