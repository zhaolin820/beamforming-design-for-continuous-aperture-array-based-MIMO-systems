function [Rate] = algorithm_WMMSE(para, H)
%The proposed WMMSE algorithm for CAPA-MIMO beamforming
%
%  [Rate] = algorithm_WMMSE(para, H)
%Inputs:
%   para: structure of the system setup
%   H: continuous CAPA-MIMO channel matrix sampled following Gauss-Legendre rule
%Outputs:
%   Rate: optimized weighted sum rate
%Date: 01/09/2025
%Author: Zhaolin Wang


%% calculate Gauss-Legendre coefficients
Phi = diag(kron(para.omega,para.omega));
Phi_T = para.Lx_T*para.Ly_T/4 * Phi; % Equation (38e)
Phi_R = para.Lx_R*para.Ly_R/4 * Phi; % Equation (38f)

%% random initialization
W = randn(para.M*para.M, para.N) + 1i*randn(para.M*para.M, para.N);
W = W./norm(W, 'fro');

%% proposed WMMSE algorithm
iter_max = 100; % maximum number of iterations
Rate_pre = 0;
for i = 1:iter_max

    sigma_2 = para.noise/para.Pt * trace(W'*Phi_T*W); % Equation (40)
    Q = W'*Phi_T*H'*Phi_R*H*Phi_T*W; % Equation (41)
    Theta = sigma_2*eye(para.N) + Q; % Theta in Equation (42)
    U = eye(para.N) + 1/sigma_2*Q; % Equation (43)
    G = Theta\W'*Phi_T*H'*Phi_R*H*Phi_T*H'*Phi_R*H*Phi_T*W/Theta; % Equation (45)
    V = Theta\W'*Phi_T*H'*Phi_R*H*Phi_T*W/Theta; % Equation (46)
    epsilon = para.Pt/para.noise/trace(U*V); % epsilon in Equation (31)
    Omega = 1/epsilon*eye(para.N) + G*U; % Omega in Equation (47)
    W = H'*Phi_R*H*Phi_T*W/Theta*U/Omega; % Equation (48)

    % calculate the rate using Equation (49)
    Rate = real(log2(det(eye(para.N) + 1/sigma_2*Q)));
    
    % check convergence
    if abs(Rate - Rate_pre)/Rate < 1e-4
        break;
    end
    Rate_pre = Rate;
end


end

