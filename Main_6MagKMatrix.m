% IMPORTANT - This script requires the Matlab symbolic toolbox
% The code is the implementations of Double-Stage Kalman.

% Author: Xinzhe Gui
% Email: lwft@qq.com
%% environment config
clear all;
reset(symengine);
load('./Mat/StatePrediction.mat');
addpath(genpath(pwd));
%% derive equations for fusion of true mag measurements
hk2 = Tnb*[1;0;0]; % predicted measurement

Hk2 = jacobian(hk2,stateVector); % measurement Jacobian
% [Hk2,SHk2]=OptimiseAlgebra(Hk2,'SHk2'); % optimise processing
SHk2=nan;

K_MagX = (Pk_k1*transpose(Hk2(1,:)))/(Hk2(1,:)*Pk_k1*transpose(Hk2(1,:)) + R_Mag);
[K_MagX,SK_MagX]=OptimiseAlgebra(K_MagX,'SK_MagX'); % Kalman gain vector
% SK_MagX=0;
K_MagY = (Pk_k1*transpose(Hk2(2,:)))/(Hk2(2,:)*Pk_k1*transpose(Hk2(2,:)) + R_Mag);
[K_MagY,SK_MagY]=OptimiseAlgebra(K_MagY,'SK_MagY'); % Kalman gain vector
% SK_MagY=0;
K_MagZ = (Pk_k1*transpose(Hk2(3,:)))/(Hk2(3,:)*Pk_k1*transpose(Hk2(3,:)) + R_Mag);
[K_MagZ,SK_MagZ]=OptimiseAlgebra(K_MagZ,'SK_MagZ'); % Kalman gain vector
% SK_MagZ=0;

% save equations and reset workspace
save('./Mat/Mag.mat','hk2','SHk2','Hk2','SK_MagX','SK_MagY','SK_MagZ','K_MagX','K_MagY','K_MagZ');