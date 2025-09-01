function [Rate] = algorithm_SPDA_SVD(para, H)
%SPDA-MIMO beamforming
%
%  [Rate] = algorithm_WMMSE(para, H)
%Inputs:
%   para: structure of the system setup
%   H: SPDA-MIMO channel matrix
%Outputs:
%   Rate: optimized weighted sum rate
%Date: 01/09/2025
%Author: Zhaolin Wang


K = H'*H; % SPDA channel kernel

% eigendecomposition of the kernel
[~, Lambda] = eig(K);
Lambda = diag(real(Lambda)); Lambda = sort(Lambda, "descend");
Lambda = abs(Lambda(1:para.N));


%% The following is the classical water filling algorithm. 
% The code refers to https://blog.csdn.net/stay_alive_13/article/details/112846063
% loss: a Kx1 vector
% P: total power to be allocated
% A: diagonal matrix, i.e. weight
% power_allocation: diagonal matrix, power allocation result, trace equals P

loss = para.noise./Lambda;
A = eye(para.N);
P = para.Pt;

K = length(loss); % number of users

w = diag(A); % width
h = loss./w; % height

allo_set = 1:K; % Initialize the index collection of users to be filled with water
level = (P+sum(loss))/sum(w); % virtual water level
[h_hat, k_hat] = max(h);

while h_hat>=level
    
    allo_set(k_hat) = -1;
    level = (P+sum(loss(allo_set>0)))/sum(w(allo_set>0)); % virtual water level
    [h_hat, ~] = max(h(allo_set>0));
    k_hat = find(h == h_hat);

end

% calculate power allocation matrix
power_allocation = zeros(K,1);
for k = 1:K
    if allo_set(k)>0
        power_allocation(k) = (level - h(k))*w(k);
    end
end


%% calculate weighted sum rate
Rate = sum(log2(1 + power_allocation./loss)); 
end

