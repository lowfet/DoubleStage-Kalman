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
 
% Derive the predicted covariance matrix using the standard equation
Pk_k1 = F*Pk*transpose(F) + Q;
% no optimise, because of the poor optimised results
% [Pk_k1,SPk_k1]=OptimiseAlgebra(Pk_k1,'SPk_k1');
SPk_k1=nan;

save('./Mat/CovariancePrediction_EasyQ.mat');