% ==============================================================================
% Copyright (c) [2024] [Zhou Bo]
% All Rights Reserved
%
% This code is the intellectual property of [Your Full Name]. The code
% or any portion of it may be used exclusively for academic and scientific
% research purposes. Any other use, including but not limited to commercial
% application or adaptation, redistribution, or incorporation into other
% software, is strictly prohibited without prior written consent from the
% copyright owner.
%
% By using this code, you agree to abide by these terms and acknowledge that
% this code is provided "AS IS" without warranty of any kind, either expressed
% or implied, including but not limited to the implied warranties of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% ==============================================================================
clear all
rng("shuffle")

NN=5000;%The number of independent repetitions of the simulation, default 5000.
R0=-1;%R_0 of the paper. Here, we provide two options, R0=-1 or R0=-0.985,...
% other inputs will be defaulted to R0=-1.
Nl=50;%N_l in the paper.
Aub=1;%A_ub in the paper, Alb is defaulted to 0.
varingCenter=1/2;%The proportion of the center position of the time-varying...
% region in the DTC. For example, 1/2 indicates that the time-varying region...
% is located exactly in the center of the DTC.
L3=1/4;%L_3 in the paper, the length of the time-varying region.L3 should...
% be a multiple of 1/128.
NP=10;%NP=T_m/T_0.
Ns=1;%N_s in the paper.
L0=1/2;%L_0 in the paper.L0 should be a multiple of 1/128.
jt=0.1;%Transition time to wave source period ratio. The transition time is...
% set to avoid numerical reflections.
rate_delta=jt/NP;
Ex=0;
if R0==-0.985
    Ex_R=zeros(3+L0*128,NN);%Ex_R is the stored E_x(z,t).
else
    Ex_R=zeros(1+L0*128,NN);%Ex_R is the stored E_x(z,t).
end
formatSpec = "%.2flam_%.2flam_varingCenter=%.2f_dfu=%.5f_NODT=%.2f_R=%.2f_NP=%.2f_TN=%.2f_NN=%.2f.mat";
        filename = compose(formatSpec,L0,...
            L3,varingCenter,Aub,Nl,R0,NP,Ns,NN);
for inx=1:NN
    [Ex,epsr_end]=FDTD_ZB_1D(Nl,...
        Aub,NP,rate_delta,Ns,L0,L3,varingCenter,R0);
    Ex_R(:,inx)=Ex;
    
    if mod(inx,100)==0
        save(filename)
        inx
    end
end
histogram(Ex_R(10,:));%Plot the histogram of the probability density of...
% E_x(z=10*lambda_0/128,t=(Ns+Nl)*T_0).
