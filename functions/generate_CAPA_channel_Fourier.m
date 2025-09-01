function [H_w] = generate_CAPA_channel_Fourier(para, H_GL)
%Calculate the CAPA-MIMO channel in the wavenumber domain
%
%  [H_w] = generate_CAPA_channel_GL(para, H)
%Inputs:
%   para: structure of the initial parameters
%   H_GL: continuous CAPA-MIMO channel matrix sampled following Gauss-Legendre rule
%Outputs:
%   H_w: CAPA-MIMO channel matrix in the wavenumber domain
%Date: 01/09/2025
%Author: Zhaolin Wang


Phi = diag(kron(para.omega,para.omega));
Phi_T = para.Lx_T*para.Ly_T/4 * Phi;
Phi_R = para.Lx_R*para.Ly_R/4 * Phi;

% number of sampling point at the transmitter
Mx_T = ceil(para.Lx_T/para.lambda); 
My_T = ceil(para.Ly_T/para.lambda);
M_T = (2*Mx_T+1)*(2*My_T+1);

% number of sampling point at the receiver
Mx_R = ceil(para.Lx_R/para.lambda); 
My_R = ceil(para.Ly_R/para.lambda);
M_R = (2*Mx_R+1)*(2*My_R+1);


%% Calculate the wavenumber-domain channel
H_w = zeros(M_R, M_T);
for p = 1:M_T
    for q = 1:M_R
        
        % wavenumber-domain sampling point at the transmitter
        n_T = ceil(p/(2*Mx_T+1));
        m_T = p - (n_T-1)*(2*My_T+1);
        n_T = n_T - (Mx_T+1);
        m_T = m_T - (My_T+1);
        Psi_T = GL_Fourier_vector(para, n_T, m_T, 'T');

        % wavenumber-domain sampling point at the receiver
        n_R = ceil(q/(2*Mx_R+1));
        m_R = q - (n_R-1)*(2*My_R+1);
        n_R = n_R - (Mx_R+1);
        m_R = m_R - (My_R+1);
        Psi_R = GL_Fourier_vector(para, n_R, m_R, 'R');

        % wavenumer-domain channel, Equation (71)
        H_w(q, p) = Psi_R'*Phi_R*H_GL*Phi_T*Psi_T;

    end
end

end


function [Psi_vec] = GL_Fourier_vector(para, n, m, mode)

Psi_vec = zeros(para.M*para.M, 1);

for p=1:para.M*para.M        
    i = ceil(p/para.M);
    j = p - (i-1)*para.M;

    if strcmp(mode, 'T')
        s_x = para.theta(i)*para.Lx_T/2;
        s_y = para.theta(j)*para.Ly_T/2;
    elseif strcmp(mode, 'R')
        s_x = para.theta(i)*para.Lx_R/2;
        s_y = para.theta(j)*para.Ly_R/2;
    end

    Psi_vec(p) = Fourier_basis_function(para, n, m, s_x, s_y, mode);
end

end