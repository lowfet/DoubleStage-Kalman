% IMPORTANT - This script requires the Matlab symbolic toolbox
% The code is the implementations of Double-Stage Kalman.

% Author: Xinzhe Gui
% Email: lwft@qq.com
%% environment config
clear all;
reset(symengine);
load('./Mat/StatePrediction.mat');
addpath(genpath(pwd));
%% derive equations for fusion of true acc measurements
hk1 = g*Tnb*[0;0;1]; % predicted measurement
% [hk1,S_hk1]=OptimiseAlgebra(hk1,'S_hk1'); % optimise processing
% S_hk1=0;
Hk1 = jacobian(hk1,stateVector); % measurement Jacobian
% [Hk1,SHk1]=OptimiseAlgebra(Hk1,'SHk1'); % optimise processing
SHk1=nan;

K_dVelX = (Pk_k1*transpose(Hk1(1,:)))/(Hk1(1,:)*Pk_k1*transpose(Hk1(1,:)) + R_dVel);
[K_dVelX,SK_dVelX]=OptimiseAlgebra(K_dVelX,'SK_dVelX'); % Kalman gain vector
% SK_dVelX=0;
K_dVelY = (Pk_k1*transpose(Hk1(2,:)))/(Hk1(2,:)*Pk_k1*transpose(Hk1(2,:)) + R_dVel);
[K_dVelY,SK_dVelY]=OptimiseAlgebra(K_dVelY,'SK_dVelY'); % Kalman gain vector
% SK_dVelY=0;
K_dVelZ = (Pk_k1*transpose(Hk1(3,:)))/(Hk1(3,:)*Pk_k1*transpose(Hk1(3,:)) + R_dVel);
[K_dVelZ,SK_dVelZ]=OptimiseAlgebra(K_dVelZ,'SK_dVelZ'); % Kalman gain vector
% SK_dVelZ=0;

% save equations and reset workspace
save('./Mat/dVel.mat','hk1','SHk1','Hk1','SK_dVelX','SK_dVelY','SK_dVelZ','K_dVelX','K_dVelY','K_dVelZ');