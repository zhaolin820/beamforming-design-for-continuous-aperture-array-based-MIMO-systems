function [H] = generate_SPDA_channel(para)
%Calculate the SPDA-MIMO channel
%
%  [H] = generate_SPDA_channel(para)
%Inputs:
%   para: structure of the initial parameters
%Outputs:
%   H: SPDA-MIMO channel matrix
%Date: 01/09/2025
%Author: Zhaolin Wang

 % half-wavelength antenna spacing
d = para.lambda/2;

% number of antennas at the transmitter
Nx_T = ceil(para.Lx_T/d);
Ny_T = ceil(para.Ly_T/d);
N_T = Nx_T*Ny_T;

% number of antennas at the receiver
Nx_R = ceil(para.Lx_R/d);
Ny_R = ceil(para.Ly_R/d);
N_R = Nx_R*Ny_R;

%% calculate channel of SPDA
H = zeros(N_R, N_T);
for p = 1:N_T
    for q = 1:N_R

        n_T = ceil(p/Nx_T);
        m_T = p - (n_T-1)*Ny_T;
        s_x = -para.Lx_T/2 + (n_T-3/2)*d;
        s_y = -para.Ly_T/2 + (m_T-3/2)*d;
        s = [s_x, s_y, 0]';


        n_R = ceil(q/Nx_R);
        m_R = q - (n_R-1)*Ny_R;
        r_x = -para.Lx_R/2 + (n_R-3/2)*d;
        r_y = -para.Ly_R/2 + (m_R-3/2)*d;
        r = para.R_phi*[r_x, r_y, 0]' + para.r_O;

        H_qp = channel_function(para, r, s);
        H(q,p)= H_qp * para.lambda^2/(4*pi); 

    end
end


end