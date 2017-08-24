% IMPORTANT - This script requires the Matlab symbolic toolbox
% The code is the implementations of Double-Stage Kalman.

% Author: Xinzhe Gui
% Email: lwft@qq.com

% State vector:
% attitude quaternion

% Time varying parameters:
% XYZ delta angle measurements in body axes - rad/s

% Observations:
% XYZ delta velocity measurements in body axes (Accerolator)- m^2/s
% Body Magnetic Field Vector - (X,Y,Z)
%% environment config
clear all;
reset(symengine);
addpath(genpath(pwd));
%% define symbolic variables and constants
syms q0 q1 q2 q3 real % quaternions defining attitude of body axes relative to local NED
syms dax day daz real % IMU delta angle measurements in body axes - rad
syms dvx dvy dvz real % IMU delta velocity measurements in body axes - m/sec
syms magX magY magZ real; % XYZ body fixed magnetic field measurements - milligauss
syms daxVar dayVar dazVar dvxVar dvyVar dvzVar real; % IMU delta angle and delta velocity measurement variances
syms dt real % IMU time step - sec
syms R_dVel R_Mag real
syms g real
%% define the state prediction equations
% define the measured Delta angle and delta velocity vectors
dAngTruth = [dax; day; daz];
dVelTruth = [dvx; dvy; dvz];

% define the quaternion rotation vector for the state estimate
quat = [q0;q1;q2;q3];
% derive the truth body to nav direction cosine matrix
Tbn = Quat2Tbn(quat);
Tnb = transpose(Tbn);

% define the attitude update equations
% use a first order expansion of rotation to calculate the quaternion increment
% acceptable for propagation of covariances
% one step state prediction
deltaQuat = [1;
    0.5*dAngTruth(1)*dt;
    0.5*dAngTruth(2)*dt;
    0.5*dAngTruth(3)*dt;
    ];
quatNew = QuatMult(quat,deltaQuat);
% the method before has the similar result with the method below
% deltaQuat1= [0;
%     dAngTruth(1);
%     dAngTruth(2);
%     dAngTruth(3);
%     ];
% quatNew1=quat+1/2*QuatMult(quat,deltaQuat1)*dt;

% Define the state vector & number of states
stateVector = quat;
nStates = numel(stateVector);
nObservers1 = 3;
nObservers2 = 3;

% Define vector of process equations
newStateVector = quatNew;

% derive the state transition matrix
F = jacobian(newStateVector, stateVector);
% the F=-F is need, because of the function(subexpr) bug in the
% function(OptimiseAlgebra)
F=-F;
[F,SF]=OptimiseAlgebra(F,'SF');
F=-F;

% % Fix The Error of Matlab symbolic toolbox
% SF =[(dax*dt)/2;
%     (daz*dt)/2;
%     (day*dt)/2];
% F =[    1, -SF(1), -SF(3), -SF(2);
%     SF(1),      1,  SF(2), -SF(3);
%     SF(3), -SF(2),      1,  SF(1);
%     SF(2),  SF(3), -SF(1),      1];

% define a symbolic covariance matrix using strings to represent 
% '_l_' to represent '( '
% '_c_' to represent ,
% '_r_' to represent ')' 
% these can be substituted later to create executable code
for rowIndex = 1:nStates
    for colIndex = 1:nStates
        eval(['syms Pk_l_',num2str(rowIndex),'_c_',num2str(colIndex), '_r_ real']);
        eval(['Pk(',num2str(rowIndex),',',num2str(colIndex), ') = Pk_l_',num2str(rowIndex),'_c_',num2str(colIndex),'_r_;']);
    end
end

for rowIndex = 1:nStates
    for colIndex = 1:nStates
        eval(['syms Pk_k1_l_',num2str(rowIndex),'_c_',num2str(colIndex), '_r_ real']);
        eval(['Pk_k1(',num2str(rowIndex),',',num2str(colIndex), ') = Pk_k1_l_',num2str(rowIndex),'_c_',num2str(colIndex),'_r_;']);
    end
end

for rowIndex = 1:nStates
    for colIndex = 1:nObservers1
        eval(['syms Kk1_l_',num2str(rowIndex),'_c_',num2str(colIndex), '_r_ real']);
        eval(['Kk1(',num2str(rowIndex),',',num2str(colIndex), ') = Kk1_l_',num2str(rowIndex),'_c_',num2str(colIndex),'_r_;']);
    end
end
for rowIndex = 1:3
    for colIndex = 1:1
        eval(['syms hk1_l_',num2str(rowIndex),'_c_',num2str(colIndex), '_r_ real']);
        eval(['hk1(',num2str(rowIndex),',',num2str(colIndex), ') = hk1_l_',num2str(rowIndex),'_c_',num2str(colIndex),'_r_;']);
    end
end
for rowIndex = 1:3
    for colIndex = 1:4
        eval(['syms Hk1_l_',num2str(rowIndex),'_c_',num2str(colIndex), '_r_ real']);
        eval(['Hk1(',num2str(rowIndex),',',num2str(colIndex), ') = Hk1_l_',num2str(rowIndex),'_c_',num2str(colIndex),'_r_;']);
    end
end
for rowIndex = 1:nStates
    for colIndex = 1:nStates
        eval(['syms Pk1_l_',num2str(rowIndex),'_c_',num2str(colIndex), '_r_ real']);
        eval(['Pk1(',num2str(rowIndex),',',num2str(colIndex), ') = Pk1_l_',num2str(rowIndex),'_c_',num2str(colIndex),'_r_;']);
    end
end

for rowIndex = 1:nStates
    for colIndex = 1:nObservers1
        eval(['syms Kk2_l_',num2str(rowIndex),'_c_',num2str(colIndex), '_r_ real']);
        eval(['Kk2(',num2str(rowIndex),',',num2str(colIndex), ') = Kk2_l_',num2str(rowIndex),'_c_',num2str(colIndex),'_r_;']);
    end
end
for rowIndex = 1:3
    for colIndex = 1:1
        eval(['syms hk2_l_',num2str(rowIndex),'_c_',num2str(colIndex), '_r_ real']);
        eval(['hk2(',num2str(rowIndex),',',num2str(colIndex), ') = hk2_l_',num2str(rowIndex),'_c_',num2str(colIndex),'_r_;']);
    end
end
for rowIndex = 1:3
    for colIndex = 1:4
        eval(['syms Hk2_l_',num2str(rowIndex),'_c_',num2str(colIndex), '_r_ real']);
        eval(['Hk2(',num2str(rowIndex),',',num2str(colIndex), ') = Hk2_l_',num2str(rowIndex),'_c_',num2str(colIndex),'_r_;']);
    end
end

%EasyQ Use
 for rowIndex = 1:4
    for colIndex = 1:4
        eval(['syms Q_l_',num2str(rowIndex),'_c_',num2str(colIndex), '_r_ real']);
        eval(['Q(',num2str(rowIndex),',',num2str(colIndex), ') = Q_l_',num2str(rowIndex),'_c_',num2str(colIndex),'_r_;']);
    end
 end

save './Mat/StatePrediction.mat';