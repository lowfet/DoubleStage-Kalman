% IMPORTANT - This script requires the Matlab symbolic toolbox
% The code is the implementations of Double-Stage Kalman.

% Author: Xinzhe Gui
% Email: lwft@qq.com
%% environment config
clear all;
reset(symengine);
load('./Mat/StatePrediction.mat');
addpath(genpath(pwd));
%% derive the covariance prediction equations
% derive the control(disturbance) influence matrix from IMu noise to state
% noise
G = jacobian(newStateVector, dAngTruth);
[G,SG]=OptimiseAlgebra(G,'SG');

% Fix The Error of Matlab symbolic toolbox
% SG =[(dt*q0)/2;
%     (dt*q1)/2;
%     (dt*q2)/2;
%     (dt*q3)/2;];
% 
% G =[ -SG(2), -SG(3), -SG(4);
%       SG(1), -SG(4),  SG(3);
%       SG(4),  SG(1), -SG(2);
%      -SG(3),  SG(2),  SG(1)];

% derive the state error matrix
distMatrix = diag([daxVar dayVar dazVar]);
Q = G*distMatrix*transpose(G);%如果Q是对角线矩阵，则就和我们本身的形式是一样的
[Q,SQ]=OptimiseAlgebra(Q,'SQ');

% SQ = dt^2;
% Q =[(daxVar*dt^2*q1^2)/4 + (dayVar*dt^2*q2^2)/4 + (dazVar*dt^2*q3^2)/4, (dayVar*dt^2*q2*q3)/4 - (daxVar*dt^2*q0*q1)/4 - (dazVar*dt^2*q2*q3)/4, (dazVar*dt^2*q1*q3)/4 - (dayVar*dt^2*q0*q2)/4 - (daxVar*dt^2*q1*q3)/4, (daxVar*dt^2*q1*q2)/4 - (dayVar*dt^2*q1*q2)/4 - (dazVar*dt^2*q0*q3)/4;
% (dayVar*dt^2*q2*q3)/4 - (daxVar*dt^2*q0*q1)/4 - (dazVar*dt^2*q2*q3)/4, (daxVar*dt^2*q0^2)/4 + (dazVar*dt^2*q2^2)/4 + (dayVar*dt^2*q3^2)/4, (daxVar*dt^2*q0*q3)/4 - (dayVar*dt^2*q0*q3)/4 - (dazVar*dt^2*q1*q2)/4, (dazVar*dt^2*q0*q2)/4 - (dayVar*dt^2*q1*q3)/4 - (daxVar*dt^2*q0*q2)/4;
% (dazVar*dt^2*q1*q3)/4 - (dayVar*dt^2*q0*q2)/4 - (daxVar*dt^2*q1*q3)/4, (daxVar*dt^2*q0*q3)/4 - (dayVar*dt^2*q0*q3)/4 - (dazVar*dt^2*q1*q2)/4,    (dayVar*dt^2*q0^2)/4 + (dazVar*dt^2*q1^2)/4 + (daxVar*dt^2*q3^2)/4, (dayVar*dt^2*q0*q1)/4 - (daxVar*dt^2*q2*q3)/4 - (dazVar*dt^2*q0*q1)/4;
% (daxVar*dt^2*q1*q2)/4 - (dayVar*dt^2*q1*q2)/4 - (dazVar*dt^2*q0*q3)/4, (dazVar*dt^2*q0*q2)/4 - (dayVar*dt^2*q1*q3)/4 - (daxVar*dt^2*q0*q2)/4, (dayVar*dt^2*q0*q1)/4 - (daxVar*dt^2*q2*q3)/4 - (dazVar*dt^2*q0*q1)/4,    (dazVar*dt^2*q0^2)/4 + (dayVar*dt^2*q1^2)/4 + (daxVar*dt^2*q2^2)/4];

% Derive the predicted covariance matrix using the standard equation
Pk_k1 = F*Pk*transpose(F) + Q;

% Collect common expressions to optimise processing
[Pk_k1,SPk_k1]=OptimiseAlgebra(Pk_k1,'SPk_k1');

save('./Mat/CovariancePrediction.mat');