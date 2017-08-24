% IMPORTANT - This script requires the Matlab symbolic toolbox
% The code is the implementations of Double-Stage Kalman.

% Author: Xinzhe Gui
% Email: lwft@qq.com
%% environment config
clear all;
reset(symengine);
load('./Mat/StatePrediction.mat');
addpath(genpath(pwd));
%% derive equations for qc2 and Pk2(Pk)
Zk2 = [magX;magY;magZ];

qc2 = Kk2*(Zk2-hk2);
[qc2,Sqc2]=OptimiseAlgebra(qc2,'Sqc2'); 
% Sqc2=0;

Pk2 = Pk1-Kk2*Hk2*Pk1;
[Pk2,SPk2]=OptimiseAlgebra(Pk2,'SPk2'); 

% save equations and reset workspace
save('./Mat/qc2AndPk2.mat','qc2','Sqc2','Pk2','SPk2');