function [gamma_e] = noise_level(STFT)
%NOISE_LEVEL Estimation of the noise level.
%   [gamma] = noise_level(STFT)
%
% INPUTS:
%   STFT    : Noisy short time fourier transform.
%
% OUTPUTS:
%   gamma_e : For a window g and a noise of variance sigma_n,
%             it is an estimation of sigma_n ||g||_2
%             see [1].
%
% REFERENCES:
% [1] D. Donoho and I. Johnstone, ???Ideal spatial adaptation via wavelet
% shrinkage,??? Biometrika, vol. 81, pp. 425???455, 1994.

    gamma_e = median(abs(real(STFT(:))))/0.6745;
end