function [p_shape, p_scale] = vor_energy_d_param()
%VOR_ENERGY_D_PARAM Parameters for the energy distribution on voronoi cells
%   Model using Gamma distribution with shape, scale parameters
%   p_shape
%   p_scale
%   below values are set using numerical estimation

p_shape = 2.284350;% k
p_scale = 1.703826;% theta
end

