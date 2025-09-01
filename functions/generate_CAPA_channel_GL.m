function [H_GL] = generate_CAPA_channel_GL(para)
%Calculate the sampled CAPA-MIMO matrix following Gauss-Legendre quadrature 
%
%  [H_GL] = generate_CAPA_channel_GL(para)
%Inputs:
%   para: structure of the initial parameters
%Outputs:
%   H_GL: continuous CAPA-MIMO channel matrix sampled following Gauss-Legendre rule
%Date: 01/09/2025
%Author: Zhaolin Wang

H_GL = zeros(para.M*para.M, para.M*para.M);
for p=1:para.M*para.M % index for Tx-CAPA
    for q=1:para.M*para.M % index for Rx-CAPA
        
        % sampling at the transmitter following Equation (34)
        i_T = ceil(p/para.M);
        j_T = p - (i_T-1)*para.M;
        s = [ para.theta(i_T)*para.Lx_T/2, para.theta(j_T)*para.Ly_T/2, 0]';

        % sampling at the receiver following Equation (37)
        i_R = ceil(q/para.M);
        j_R = q - (i_R-1)*para.M;
        r = para.R_phi*[ para.theta(i_R)*para.Lx_R/2, para.theta(j_R)*para.Ly_R/2, 0]' + para.r_O;
        
        % calculate the channel between the sampled points
        H_GL(q,p) = channel_function(para, r, s);
    end
end

end