function [H] = channel_function(para, r, s)
%Calculate the channel between the point s on Tx-CAPA and the point r on Rx-CAPA
%
%  [H] = channel_function(para, r, s)
%Inputs:
%   para: structure of the system setup
%   r: a point on Rx-CAPA
%   s: a point on Tx-CAPA
%Outputs:
%   H: channel
%Date: 01/09/2025
%Author: Zhaolin Wang

d = norm(r-s); % distance between the two points

u_T = [0,1,0]'; % transmit polarization direction
u_R = para.R_phi*u_T; % receive polarization direction

H = -1i*para.eta*exp(-1i*2*pi/para.lambda*d)/(2*para.lambda*d)...
    *(eye(3) - (r-s)*(r-s)'/d^2); % Green's function in Equation (3)
H = u_R'*H*u_T; % uni-polarized channel



end

