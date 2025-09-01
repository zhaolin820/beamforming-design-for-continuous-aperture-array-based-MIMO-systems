function [Phi] = Fourier_basis_function(para, n_x, n_y, s_x, s_y, mode)
%Fourier basis function for a point (s_x, s_y, 0) on CAPA
%
%  [Phi] = Fourier_basis_function(para, n_x, n_y, s_x, s_y)
%Inputs:
%   para: structure of the initial parameters
%   [n_x, n_y]: index of the Fourier basis function
%   [s_x, s_y]: x- and y-coordinates of a point on CAPA
%Outputs:
%   Phi: Fourier basis functions
%Date: 06/03/2025
%Author: Zhaolin Wang

if strcmp(mode, 'T')
    Lx = para.Lx_T;
    Ly = para.Ly_T;
elseif strcmp(mode, 'R')
    Lx = para.Lx_R;
    Ly = para.Ly_R;
end

A = Lx*Ly;
Phi = 1/sqrt(A) * exp(1i*2*pi * ( n_x/Lx.*s_x + n_y/Ly.*s_y ) );

end

