% IMPORTANT - This script requires the Matlab symbolic toolbox
% The code is the implementations of Double-Stage Kalman.

% Author: Xinzhe Gui
% Email: lwft@qq.com

%% environment config
clear all;
reset(symengine);
load('./Mat/StatePrediction.mat');
addpath(genpath(pwd));
%% derive equations for qc1 and Pk1
Zk1 = [dvx;dvy;dvz];

qc1 = Kk1*(Zk1-hk1);
[qc1,Sqc1]=OptimiseAlgebra(qc1,'Sqc1');
% Sqc1=0;

Pk1 = Pk_k1-Kk1*Hk1*Pk_k1;
[Pk1,SPk1]=OptimiseAlgebra(Pk1,'SPk1');

% save equations and reset workspace
save('./Mat/qc1AndPk1.mat','qc1','Sqc1','Pk1','SPk1');