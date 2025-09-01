clc
clear all
close all

addpath("functions");

%% System setup
para = para_init();

%% CAPA, Proposed algorithm
[H] = generate_CAPA_channel_GL(para);
rate_CAPA = algorithm_WMMSE(para, H);
disp(['CAPA, Proposed WMMSE: Rate - ' num2str(rate_CAPA) ' bit/s/Hz']);

%% CAPA, Fourier-based approach
H_w = generate_CAPA_channel_Fourier(para, H);
rate_CAPA_Fourier = algorithm_Fourier_SVD(para, H_w);
disp(['CAPA, Fourier-SVD: Rate - ' num2str(rate_CAPA_Fourier) ' bit/s/Hz']);

%% SPDA
H_SPDA = generate_SPDA_channel(para);
rate_SPDA = algorithm_SPDA_SVD(para, H_SPDA);
disp(['SPDA, SVD: Rate - ' num2str(rate_SPDA) ' bit/s/Hz']);
