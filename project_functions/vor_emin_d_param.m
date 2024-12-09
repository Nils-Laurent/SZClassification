function [p_shape, p_scale] = vor_emin_d_param()
%VOR_ENERGY_D_PARAM Parameters for the energy distribution on voronoi cells
%   Model using Gamma distribution with shape, scale parameters
%   p_shape
%   p_scale
%   below values are set using numerical estimation

p_shape = 1.12623;% k
p_scale = 1.83190;% theta
end

