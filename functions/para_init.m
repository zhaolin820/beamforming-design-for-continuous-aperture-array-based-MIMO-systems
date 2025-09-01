function [para] = para_init()
%Construct a struct of the system parameters 
%  [para] = para_init()
%Inputs:
%   None
%Outputs:
%   para: a struct
%Date: 01/09/2025
%Author: Zhaolin Wang

para.fc = 2.4e9; % carrier frequency
c = 3e8; % speed of light
para.lambda = c/para.fc; % wavelength

para.Pt = 0.1; % overall transmit power (A^2)
para.noise = 5.6e-3; % noise power (V^2/m^2)

para.N = 10; % number of data streams
para.M = 20; % number of GL samples

% aperture area
para.A_T = 0.25; % m^2
para.A_R = 0.25; % m^2

% aperture region
para.Lx_T = sqrt(para.A_T);
para.Ly_T = sqrt(para.A_T);
para.Lx_R = sqrt(para.A_R);
para.Ly_R = sqrt(para.A_R);

% location of the receiver
para.r_O = [0,0,10]';

% rotation of the receiver
para.phi = 0/180*pi;
para.R_phi = [1, 0, 0; 
    0, cos(para.phi), -sin(para.phi); 
    0, sin(para.phi), cos(para.phi)];

% free-space impedance
para.eta = 120*pi;

% Gauss-Legendre quadrature parameters
[theta, omega] = GaussLegendre(para.M);
para.theta = theta;
para.omega = omega;

end

