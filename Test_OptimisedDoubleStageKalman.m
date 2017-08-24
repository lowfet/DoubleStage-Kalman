% IMPORTANT - This script requires the Matlab symbolic toolbox
% The code is the implementations of Double-Stage Kalman.

% Author: Xinzhe Gui, Jing Xue
% Email: lwft@qq.com
%%
clear all;
close all;

addpath(genpath(pwd));

AccSensx=0.995544148719357;
AccSensy=1.00067100989279;
AccSensz=0.999918078434263;
AccOffsetx=0.00428358112745599;
AccOffsety=-0.00270633319254957;
AccOffsetz=-0.0705557225725813;

magcalibration=[2.06385073217840,0.0413407668037351,0.0246229991758043;0,2.10288949976811,-0.0104899029262219;0,0,2.36140472565193];
magoffset=[0.0610603654152471;0.0201527841522371;0.990707291320690];

Gyro_x_offset=-19.602; Gyro_y_offset=11.886; Gyro_z_offset=1.338;

IMU=xlsread('TestData/datal_1.xlsx','Sheet2');
[mIMUData,nIMUData]=size(IMU);
N=mIMUData;

g=9.8;unit=32768;Mag_unit=0.3;
q=[1;0;0;0];
% q=angle2quat(IMU(1,12)*pi/180,IMU(1,11)*pi/180,IMU(1,10)*pi/180)';
Pk=eye(4,4)*(0.125-0.0003)+ones(4,4)*0.0003;
daxVar=0.000001;dayVar=0.000001;dazVar=0.000001;
R_dVel=2;
R_Mag=1;
Q=10^(-6)*eye(4,4);
time=0;
for k=1:N
    
    normq=sqrt(q(1)^2+q(2)^2+q(3)^2+q(4)^2);
    q(1)=q(1)/normq;
    q(2)=q(2)/normq;
    q(3)=q(3)/normq;
    q(4)=q(4)/normq;
    
    Cnb = [q(1)^2+q(2)^2-(q(3)^2+q(4)^2),     2*(q(2)*q(3)+q(1)*q(4)),     2*(q(2)*q(4)-q(1)*q(3));
       2*(q(2)*q(3)-q(1)*q(4)),         q(1)^2+q(3)^2-(q(2)^2+q(4)^2),     2*(q(1)*q(2)+q(3)*q(4));
       2*(q(2)*q(4)+q(1)*q(3)),     2*(q(3)*q(4)-q(1)*q(2)),        q(1)^2+q(4)^2-(q(2)^2+q(3)^2)];

    Phi(k)=atan2(2*(q(3)*q(4)+q(1)*q(2)),q(1)^2+q(4)^2-(q(2)^2+q(3)^2))/pi*180;%ROLL*pi/180  
    Theta(k)=asin(2*(q(1)*q(3)-q(4)*q(2)))/pi*180;%PITCH
    Psi(k)= atan2(2*(q(2)*q(3)+q(1)*q(4)),(q(1)^2+q(2)^2-q(3)^2-q(4)^2))/pi*180;%YAW

    tic
    
    wx=(IMU(k,4)-Gyro_x_offset)*2000/unit*pi/180;
    wy=-(IMU(k,5)-Gyro_y_offset)*2000/unit*pi/180;
    wz=-(IMU(k,6)-Gyro_z_offset)*2000/unit*pi/180;
    deltet=IMU(k,17)/1000000;
    
    amx=-((IMU(k,1)*8/unit)-AccOffsetx)*AccSensx*g;
    amy=((IMU(k,2)*8/unit)-AccOffsety)*AccSensy*g;
    amz=((IMU(k,3)*8/unit)-AccOffsetz)*AccSensz*g;
    
    tmpmagfx = (IMU(k,7)/100*Mag_unit);
    tmpmagfy = (IMU(k,8)/100*Mag_unit);
    tmpmagfz = (IMU(k,9)/100*Mag_unit);
    tmpmag=magcalibration*([tmpmagfx;tmpmagfy;tmpmagfz]-magoffset);
    mx = tmpmag(1);
    my = -tmpmag(2);
    mz = -tmpmag(3);
    square_sum=sqrt(mx^2+my^2+mz^2);
    mx_s=mx/square_sum;
    my_s=my/square_sum;
    mz_s=mz/square_sum;
    
    dax=wx;day=wy;daz=wz;dt=deltet;
    dvx=amx;dvy=amy;dvz=amz;
    magX=mx_s;magY=my_s;magZ=mz_s;
    
    q0=q(1);q1=q(2);q2=q(3);q3=q(4);
    
    q(1) = q0 - (dax*dt*q1)/2 - (day*dt*q2)/2 - (daz*dt*q3)/2;
    q(2) = q1 + (dax*dt*q0)/2 - (day*dt*q3)/2 + (daz*dt*q2)/2;
    q(3) = q2 + (dax*dt*q3)/2 + (day*dt*q0)/2 - (daz*dt*q1)/2;
    q(4) = q3 - (dax*dt*q2)/2 + (day*dt*q1)/2 + (daz*dt*q0)/2;
    
    q0=q(1);q1=q(2);q2=q(3);q3=q(4);
    
    SF = zeros(3,1);
    SF(1) = (day*dt)/2;
    SF(2) = (daz*dt)/2;
    SF(3) = (dax*dt)/2;
    Pk_k1 = zeros(4,4);
    Pk_k1(1,1) = Pk(1,1) + Q(1,1) - Pk(2,1)*SF(3) - Pk(3,1)*SF(1) - Pk(4,1)*SF(2) + SF(3)*(Pk(2,2)*SF(3) - Pk(1,2) + Pk(3,2)*SF(1) + Pk(4,2)*SF(2)) + SF(1)*(Pk(2,3)*SF(3) - Pk(1,3) + Pk(3,3)*SF(1) + Pk(4,3)*SF(2)) + SF(2)*(Pk(2,4)*SF(3) - Pk(1,4) + Pk(3,4)*SF(1) + Pk(4,4)*SF(2));
    Pk_k1(1,2) = Pk(1,2) + Q(1,2) - Pk(2,2)*SF(3) - Pk(3,2)*SF(1) - Pk(4,2)*SF(2) - SF(3)*(Pk(2,1)*SF(3) - Pk(1,1) + Pk(3,1)*SF(1) + Pk(4,1)*SF(2)) - SF(2)*(Pk(2,3)*SF(3) - Pk(1,3) + Pk(3,3)*SF(1) + Pk(4,3)*SF(2)) + SF(1)*(Pk(2,4)*SF(3) - Pk(1,4) + Pk(3,4)*SF(1) + Pk(4,4)*SF(2));
    Pk_k1(1,3) = Pk(1,3) + Q(1,3) - Pk(2,3)*SF(3) - Pk(3,3)*SF(1) - Pk(4,3)*SF(2) - SF(1)*(Pk(2,1)*SF(3) - Pk(1,1) + Pk(3,1)*SF(1) + Pk(4,1)*SF(2)) + SF(2)*(Pk(2,2)*SF(3) - Pk(1,2) + Pk(3,2)*SF(1) + Pk(4,2)*SF(2)) - SF(3)*(Pk(2,4)*SF(3) - Pk(1,4) + Pk(3,4)*SF(1) + Pk(4,4)*SF(2));
    Pk_k1(1,4) = Pk(1,4) + Q(1,4) - Pk(2,4)*SF(3) - Pk(3,4)*SF(1) - Pk(4,4)*SF(2) - SF(2)*(Pk(2,1)*SF(3) - Pk(1,1) + Pk(3,1)*SF(1) + Pk(4,1)*SF(2)) - SF(1)*(Pk(2,2)*SF(3) - Pk(1,2) + Pk(3,2)*SF(1) + Pk(4,2)*SF(2)) + SF(3)*(Pk(2,3)*SF(3) - Pk(1,3) + Pk(3,3)*SF(1) + Pk(4,3)*SF(2));
    Pk_k1(2,1) = Pk(2,1) + Q(2,1) + Pk(1,1)*SF(3) + Pk(3,1)*SF(2) - Pk(4,1)*SF(1) - SF(3)*(Pk(2,2) + Pk(1,2)*SF(3) + Pk(3,2)*SF(2) - Pk(4,2)*SF(1)) - SF(1)*(Pk(2,3) + Pk(1,3)*SF(3) + Pk(3,3)*SF(2) - Pk(4,3)*SF(1)) - SF(2)*(Pk(2,4) + Pk(1,4)*SF(3) + Pk(3,4)*SF(2) - Pk(4,4)*SF(1));
    Pk_k1(2,2) = Pk(2,2) + Q(2,2) + Pk(1,2)*SF(3) + Pk(3,2)*SF(2) - Pk(4,2)*SF(1) + SF(3)*(Pk(2,1) + Pk(1,1)*SF(3) + Pk(3,1)*SF(2) - Pk(4,1)*SF(1)) + SF(2)*(Pk(2,3) + Pk(1,3)*SF(3) + Pk(3,3)*SF(2) - Pk(4,3)*SF(1)) - SF(1)*(Pk(2,4) + Pk(1,4)*SF(3) + Pk(3,4)*SF(2) - Pk(4,4)*SF(1));
    Pk_k1(2,3) = Pk(2,3) + Q(2,3) + Pk(1,3)*SF(3) + Pk(3,3)*SF(2) - Pk(4,3)*SF(1) + SF(1)*(Pk(2,1) + Pk(1,1)*SF(3) + Pk(3,1)*SF(2) - Pk(4,1)*SF(1)) - SF(2)*(Pk(2,2) + Pk(1,2)*SF(3) + Pk(3,2)*SF(2) - Pk(4,2)*SF(1)) + SF(3)*(Pk(2,4) + Pk(1,4)*SF(3) + Pk(3,4)*SF(2) - Pk(4,4)*SF(1));
    Pk_k1(2,4) = Pk(2,4) + Q(2,4) + Pk(1,4)*SF(3) + Pk(3,4)*SF(2) - Pk(4,4)*SF(1) + SF(2)*(Pk(2,1) + Pk(1,1)*SF(3) + Pk(3,1)*SF(2) - Pk(4,1)*SF(1)) + SF(1)*(Pk(2,2) + Pk(1,2)*SF(3) + Pk(3,2)*SF(2) - Pk(4,2)*SF(1)) - SF(3)*(Pk(2,3) + Pk(1,3)*SF(3) + Pk(3,3)*SF(2) - Pk(4,3)*SF(1));
    Pk_k1(3,1) = Pk(3,1) + Q(3,1) + Pk(1,1)*SF(1) - Pk(2,1)*SF(2) + Pk(4,1)*SF(3) - SF(3)*(Pk(3,2) + Pk(1,2)*SF(1) - Pk(2,2)*SF(2) + Pk(4,2)*SF(3)) - SF(1)*(Pk(3,3) + Pk(1,3)*SF(1) - Pk(2,3)*SF(2) + Pk(4,3)*SF(3)) - SF(2)*(Pk(3,4) + Pk(1,4)*SF(1) - Pk(2,4)*SF(2) + Pk(4,4)*SF(3));
    Pk_k1(3,2) = Pk(3,2) + Q(3,2) + Pk(1,2)*SF(1) - Pk(2,2)*SF(2) + Pk(4,2)*SF(3) + SF(3)*(Pk(3,1) + Pk(1,1)*SF(1) - Pk(2,1)*SF(2) + Pk(4,1)*SF(3)) + SF(2)*(Pk(3,3) + Pk(1,3)*SF(1) - Pk(2,3)*SF(2) + Pk(4,3)*SF(3)) - SF(1)*(Pk(3,4) + Pk(1,4)*SF(1) - Pk(2,4)*SF(2) + Pk(4,4)*SF(3));
    Pk_k1(3,3) = Pk(3,3) + Q(3,3) + Pk(1,3)*SF(1) - Pk(2,3)*SF(2) + Pk(4,3)*SF(3) + SF(1)*(Pk(3,1) + Pk(1,1)*SF(1) - Pk(2,1)*SF(2) + Pk(4,1)*SF(3)) - SF(2)*(Pk(3,2) + Pk(1,2)*SF(1) - Pk(2,2)*SF(2) + Pk(4,2)*SF(3)) + SF(3)*(Pk(3,4) + Pk(1,4)*SF(1) - Pk(2,4)*SF(2) + Pk(4,4)*SF(3));
    Pk_k1(3,4) = Pk(3,4) + Q(3,4) + Pk(1,4)*SF(1) - Pk(2,4)*SF(2) + Pk(4,4)*SF(3) + SF(2)*(Pk(3,1) + Pk(1,1)*SF(1) - Pk(2,1)*SF(2) + Pk(4,1)*SF(3)) + SF(1)*(Pk(3,2) + Pk(1,2)*SF(1) - Pk(2,2)*SF(2) + Pk(4,2)*SF(3)) - SF(3)*(Pk(3,3) + Pk(1,3)*SF(1) - Pk(2,3)*SF(2) + Pk(4,3)*SF(3));
    Pk_k1(4,1) = Pk(4,1) + Q(4,1) + Pk(1,1)*SF(2) + Pk(2,1)*SF(1) - Pk(3,1)*SF(3) - SF(3)*(Pk(4,2) + Pk(1,2)*SF(2) + Pk(2,2)*SF(1) - Pk(3,2)*SF(3)) - SF(1)*(Pk(4,3) + Pk(1,3)*SF(2) + Pk(2,3)*SF(1) - Pk(3,3)*SF(3)) - SF(2)*(Pk(4,4) + Pk(1,4)*SF(2) + Pk(2,4)*SF(1) - Pk(3,4)*SF(3));
    Pk_k1(4,2) = Pk(4,2) + Q(4,2) + Pk(1,2)*SF(2) + Pk(2,2)*SF(1) - Pk(3,2)*SF(3) + SF(3)*(Pk(4,1) + Pk(1,1)*SF(2) + Pk(2,1)*SF(1) - Pk(3,1)*SF(3)) + SF(2)*(Pk(4,3) + Pk(1,3)*SF(2) + Pk(2,3)*SF(1) - Pk(3,3)*SF(3)) - SF(1)*(Pk(4,4) + Pk(1,4)*SF(2) + Pk(2,4)*SF(1) - Pk(3,4)*SF(3));
    Pk_k1(4,3) = Pk(4,3) + Q(4,3) + Pk(1,3)*SF(2) + Pk(2,3)*SF(1) - Pk(3,3)*SF(3) + SF(1)*(Pk(4,1) + Pk(1,1)*SF(2) + Pk(2,1)*SF(1) - Pk(3,1)*SF(3)) - SF(2)*(Pk(4,2) + Pk(1,2)*SF(2) + Pk(2,2)*SF(1) - Pk(3,2)*SF(3)) + SF(3)*(Pk(4,4) + Pk(1,4)*SF(2) + Pk(2,4)*SF(1) - Pk(3,4)*SF(3));
    Pk_k1(4,4) = Pk(4,4) + Q(4,4) + Pk(1,4)*SF(2) + Pk(2,4)*SF(1) - Pk(3,4)*SF(3) + SF(2)*(Pk(4,1) + Pk(1,1)*SF(2) + Pk(2,1)*SF(1) - Pk(3,1)*SF(3)) + SF(1)*(Pk(4,2) + Pk(1,2)*SF(2) + Pk(2,2)*SF(1) - Pk(3,2)*SF(3)) - SF(3)*(Pk(4,3) + Pk(1,3)*SF(2) + Pk(2,3)*SF(1) - Pk(3,3)*SF(3));
    Hk1 = zeros(3,4);
    Hk1(1,1) = -2*g*q2;
    Hk1(1,2) = 2*g*q3;
    Hk1(1,3) = -2*g*q0;
    Hk1(1,4) = 2*g*q1;
    Hk1(2,1) = 2*g*q1;
    Hk1(2,2) = 2*g*q0;
    Hk1(2,3) = 2*g*q3;
    Hk1(2,4) = 2*g*q2;
    Hk1(3,1) = 2*g*q0;
    Hk1(3,2) = -2*g*q1;
    Hk1(3,3) = -2*g*q2;
    Hk1(3,4) = 2*g*q3;
    SK_dVelX = zeros(1,1);
    SK_dVelX(1) = 1/(R_dVel + 2*g*q2*(2*Pk_k1(1,1)*g*q2 - 2*Pk_k1(2,1)*g*q3 + 2*Pk_k1(3,1)*g*q0 - 2*Pk_k1(4,1)*g*q1) - 2*g*q3*(2*Pk_k1(1,2)*g*q2 - 2*Pk_k1(2,2)*g*q3 + 2*Pk_k1(3,2)*g*q0 - 2*Pk_k1(4,2)*g*q1) + 2*g*q0*(2*Pk_k1(1,3)*g*q2 - 2*Pk_k1(2,3)*g*q3 + 2*Pk_k1(3,3)*g*q0 - 2*Pk_k1(4,3)*g*q1) - 2*g*q1*(2*Pk_k1(1,4)*g*q2 - 2*Pk_k1(2,4)*g*q3 + 2*Pk_k1(3,4)*g*q0 - 2*Pk_k1(4,4)*g*q1));
    SK_dVelY = zeros(1,1);
    SK_dVelY(1) = 1/(R_dVel + 2*g*q1*(2*Pk_k1(1,1)*g*q1 + 2*Pk_k1(2,1)*g*q0 + 2*Pk_k1(3,1)*g*q3 + 2*Pk_k1(4,1)*g*q2) + 2*g*q0*(2*Pk_k1(1,2)*g*q1 + 2*Pk_k1(2,2)*g*q0 + 2*Pk_k1(3,2)*g*q3 + 2*Pk_k1(4,2)*g*q2) + 2*g*q3*(2*Pk_k1(1,3)*g*q1 + 2*Pk_k1(2,3)*g*q0 + 2*Pk_k1(3,3)*g*q3 + 2*Pk_k1(4,3)*g*q2) + 2*g*q2*(2*Pk_k1(1,4)*g*q1 + 2*Pk_k1(2,4)*g*q0 + 2*Pk_k1(3,4)*g*q3 + 2*Pk_k1(4,4)*g*q2));
    SK_dVelZ = zeros(1,1);
    SK_dVelZ(1) = 1/(R_dVel + 2*g*q0*(2*Pk_k1(1,1)*g*q0 - 2*Pk_k1(2,1)*g*q1 - 2*Pk_k1(3,1)*g*q2 + 2*Pk_k1(4,1)*g*q3) - 2*g*q1*(2*Pk_k1(1,2)*g*q0 - 2*Pk_k1(2,2)*g*q1 - 2*Pk_k1(3,2)*g*q2 + 2*Pk_k1(4,2)*g*q3) - 2*g*q2*(2*Pk_k1(1,3)*g*q0 - 2*Pk_k1(2,3)*g*q1 - 2*Pk_k1(3,3)*g*q2 + 2*Pk_k1(4,3)*g*q3) + 2*g*q3*(2*Pk_k1(1,4)*g*q0 - 2*Pk_k1(2,4)*g*q1 - 2*Pk_k1(3,4)*g*q2 + 2*Pk_k1(4,4)*g*q3));
    Kk1 = zeros(4,3);
    Kk1(1,1) = -SK_dVelX(1)*(2*Pk_k1(1,1)*g*q2 + 2*Pk_k1(1,3)*g*q0 - 2*Pk_k1(1,2)*g*q3 - 2*Pk_k1(1,4)*g*q1);
    Kk1(1,2) = SK_dVelY(1)*(2*Pk_k1(1,1)*g*q1 + 2*Pk_k1(1,2)*g*q0 + 2*Pk_k1(1,3)*g*q3 + 2*Pk_k1(1,4)*g*q2);
    Kk1(1,3) = SK_dVelZ(1)*(2*Pk_k1(1,1)*g*q0 - 2*Pk_k1(1,2)*g*q1 - 2*Pk_k1(1,3)*g*q2 + 2*Pk_k1(1,4)*g*q3);
    Kk1(2,1) = -SK_dVelX(1)*(2*Pk_k1(2,1)*g*q2 + 2*Pk_k1(2,3)*g*q0 - 2*Pk_k1(2,2)*g*q3 - 2*Pk_k1(2,4)*g*q1);
    Kk1(2,2) = SK_dVelY(1)*(2*Pk_k1(2,1)*g*q1 + 2*Pk_k1(2,2)*g*q0 + 2*Pk_k1(2,3)*g*q3 + 2*Pk_k1(2,4)*g*q2);
    Kk1(2,3) = SK_dVelZ(1)*(2*Pk_k1(2,1)*g*q0 - 2*Pk_k1(2,2)*g*q1 - 2*Pk_k1(2,3)*g*q2 + 2*Pk_k1(2,4)*g*q3);
    Kk1(3,1) = -SK_dVelX(1)*(2*Pk_k1(3,1)*g*q2 + 2*Pk_k1(3,3)*g*q0 - 2*Pk_k1(3,2)*g*q3 - 2*Pk_k1(3,4)*g*q1);
    Kk1(3,2) = SK_dVelY(1)*(2*Pk_k1(3,1)*g*q1 + 2*Pk_k1(3,2)*g*q0 + 2*Pk_k1(3,3)*g*q3 + 2*Pk_k1(3,4)*g*q2);
    Kk1(3,3) = SK_dVelZ(1)*(2*Pk_k1(3,1)*g*q0 - 2*Pk_k1(3,2)*g*q1 - 2*Pk_k1(3,3)*g*q2 + 2*Pk_k1(3,4)*g*q3);
    Kk1(4,1) = -SK_dVelX(1)*(2*Pk_k1(4,1)*g*q2 + 2*Pk_k1(4,3)*g*q0 - 2*Pk_k1(4,2)*g*q3 - 2*Pk_k1(4,4)*g*q1);
    Kk1(4,2) = SK_dVelY(1)*(2*Pk_k1(4,1)*g*q1 + 2*Pk_k1(4,2)*g*q0 + 2*Pk_k1(4,3)*g*q3 + 2*Pk_k1(4,4)*g*q2);
    Kk1(4,3) = SK_dVelZ(1)*(2*Pk_k1(4,1)*g*q0 - 2*Pk_k1(4,2)*g*q1 - 2*Pk_k1(4,3)*g*q2 + 2*Pk_k1(4,4)*g*q3);
    hk1 = zeros(3,1);
    hk1(1) = -g*(2*q0*q2 - 2*q1*q3);
    hk1(2) = g*(2*q0*q1 + 2*q2*q3);
    hk1(3) = g*(q0^2 - q1^2 - q2^2 + q3^2);
    Sqc1 = zeros(3,1);
    Sqc1(1) = dvz - hk1(3,1);
    Sqc1(2) = dvy - hk1(2,1);
    Sqc1(3) = dvx - hk1(1,1);
    qc1 = zeros(4,1);
    qc1(1) = Kk1(1,1)*Sqc1(3) + Kk1(1,2)*Sqc1(2) + Kk1(1,3)*Sqc1(1);
    qc1(2) = Kk1(2,1)*Sqc1(3) + Kk1(2,2)*Sqc1(2) + Kk1(2,3)*Sqc1(1);
    qc1(3) = Kk1(3,1)*Sqc1(3) + Kk1(3,2)*Sqc1(2) + Kk1(3,3)*Sqc1(1);
    qc1(4) = Kk1(4,1)*Sqc1(3) + Kk1(4,2)*Sqc1(2) + Kk1(4,3)*Sqc1(1);
    SPk1 = zeros(16,1);
    SPk1(1) = Hk1(1,4)*Kk1(4,1) + Hk1(2,4)*Kk1(4,2) + Hk1(3,4)*Kk1(4,3);
    SPk1(2) = Hk1(1,3)*Kk1(4,1) + Hk1(2,3)*Kk1(4,2) + Hk1(3,3)*Kk1(4,3);
    SPk1(3) = Hk1(1,2)*Kk1(4,1) + Hk1(2,2)*Kk1(4,2) + Hk1(3,2)*Kk1(4,3);
    SPk1(4) = Hk1(1,1)*Kk1(4,1) + Hk1(2,1)*Kk1(4,2) + Hk1(3,1)*Kk1(4,3);
    SPk1(5) = Hk1(1,4)*Kk1(3,1) + Hk1(2,4)*Kk1(3,2) + Hk1(3,4)*Kk1(3,3);
    SPk1(6) = Hk1(1,3)*Kk1(3,1) + Hk1(2,3)*Kk1(3,2) + Hk1(3,3)*Kk1(3,3);
    SPk1(7) = Hk1(1,2)*Kk1(3,1) + Hk1(2,2)*Kk1(3,2) + Hk1(3,2)*Kk1(3,3);
    SPk1(8) = Hk1(1,1)*Kk1(3,1) + Hk1(2,1)*Kk1(3,2) + Hk1(3,1)*Kk1(3,3);
    SPk1(9) = Hk1(1,4)*Kk1(2,1) + Hk1(2,4)*Kk1(2,2) + Hk1(3,4)*Kk1(2,3);
    SPk1(10) = Hk1(1,3)*Kk1(2,1) + Hk1(2,3)*Kk1(2,2) + Hk1(3,3)*Kk1(2,3);
    SPk1(11) = Hk1(1,2)*Kk1(2,1) + Hk1(2,2)*Kk1(2,2) + Hk1(3,2)*Kk1(2,3);
    SPk1(12) = Hk1(1,1)*Kk1(2,1) + Hk1(2,1)*Kk1(2,2) + Hk1(3,1)*Kk1(2,3);
    SPk1(13) = Hk1(1,4)*Kk1(1,1) + Hk1(2,4)*Kk1(1,2) + Hk1(3,4)*Kk1(1,3);
    SPk1(14) = Hk1(1,3)*Kk1(1,1) + Hk1(2,3)*Kk1(1,2) + Hk1(3,3)*Kk1(1,3);
    SPk1(15) = Hk1(1,2)*Kk1(1,1) + Hk1(2,2)*Kk1(1,2) + Hk1(3,2)*Kk1(1,3);
    SPk1(16) = Hk1(1,1)*Kk1(1,1) + Hk1(2,1)*Kk1(1,2) + Hk1(3,1)*Kk1(1,3);
    Pk1 = zeros(4,4);
    Pk1(1,1) = Pk_k1(1,1) - Pk_k1(1,1)*SPk1(16) - Pk_k1(2,1)*SPk1(15) - Pk_k1(3,1)*SPk1(14) - Pk_k1(4,1)*SPk1(13);
    Pk1(1,2) = Pk_k1(1,2) - Pk_k1(1,2)*SPk1(16) - Pk_k1(2,2)*SPk1(15) - Pk_k1(3,2)*SPk1(14) - Pk_k1(4,2)*SPk1(13);
    Pk1(1,3) = Pk_k1(1,3) - Pk_k1(1,3)*SPk1(16) - Pk_k1(2,3)*SPk1(15) - Pk_k1(3,3)*SPk1(14) - Pk_k1(4,3)*SPk1(13);
    Pk1(1,4) = Pk_k1(1,4) - Pk_k1(1,4)*SPk1(16) - Pk_k1(2,4)*SPk1(15) - Pk_k1(3,4)*SPk1(14) - Pk_k1(4,4)*SPk1(13);
    Pk1(2,1) = Pk_k1(2,1) - Pk_k1(1,1)*SPk1(12) - Pk_k1(2,1)*SPk1(11) - Pk_k1(3,1)*SPk1(10) - Pk_k1(4,1)*SPk1(9);
    Pk1(2,2) = Pk_k1(2,2) - Pk_k1(1,2)*SPk1(12) - Pk_k1(2,2)*SPk1(11) - Pk_k1(3,2)*SPk1(10) - Pk_k1(4,2)*SPk1(9);
    Pk1(2,3) = Pk_k1(2,3) - Pk_k1(1,3)*SPk1(12) - Pk_k1(2,3)*SPk1(11) - Pk_k1(3,3)*SPk1(10) - Pk_k1(4,3)*SPk1(9);
    Pk1(2,4) = Pk_k1(2,4) - Pk_k1(1,4)*SPk1(12) - Pk_k1(2,4)*SPk1(11) - Pk_k1(3,4)*SPk1(10) - Pk_k1(4,4)*SPk1(9);
    Pk1(3,1) = Pk_k1(3,1) - Pk_k1(1,1)*SPk1(8) - Pk_k1(2,1)*SPk1(7) - Pk_k1(3,1)*SPk1(6) - Pk_k1(4,1)*SPk1(5);
    Pk1(3,2) = Pk_k1(3,2) - Pk_k1(1,2)*SPk1(8) - Pk_k1(2,2)*SPk1(7) - Pk_k1(3,2)*SPk1(6) - Pk_k1(4,2)*SPk1(5);
    Pk1(3,3) = Pk_k1(3,3) - Pk_k1(1,3)*SPk1(8) - Pk_k1(2,3)*SPk1(7) - Pk_k1(3,3)*SPk1(6) - Pk_k1(4,3)*SPk1(5);
    Pk1(3,4) = Pk_k1(3,4) - Pk_k1(1,4)*SPk1(8) - Pk_k1(2,4)*SPk1(7) - Pk_k1(3,4)*SPk1(6) - Pk_k1(4,4)*SPk1(5);
    Pk1(4,1) = Pk_k1(4,1) - Pk_k1(1,1)*SPk1(4) - Pk_k1(2,1)*SPk1(3) - Pk_k1(3,1)*SPk1(2) - Pk_k1(4,1)*SPk1(1);
    Pk1(4,2) = Pk_k1(4,2) - Pk_k1(1,2)*SPk1(4) - Pk_k1(2,2)*SPk1(3) - Pk_k1(3,2)*SPk1(2) - Pk_k1(4,2)*SPk1(1);
    Pk1(4,3) = Pk_k1(4,3) - Pk_k1(1,3)*SPk1(4) - Pk_k1(2,3)*SPk1(3) - Pk_k1(3,3)*SPk1(2) - Pk_k1(4,3)*SPk1(1);
    Pk1(4,4) = Pk_k1(4,4) - Pk_k1(1,4)*SPk1(4) - Pk_k1(2,4)*SPk1(3) - Pk_k1(3,4)*SPk1(2) - Pk_k1(4,4)*SPk1(1);
    Hk2 = zeros(3,4);
    Hk2(1,1) = 2*q0;
    Hk2(1,2) = 2*q1;
    Hk2(1,3) = -2*q2;
    Hk2(1,4) = -2*q3;
    Hk2(2,1) = -2*q3;
    Hk2(2,2) = 2*q2;
    Hk2(2,3) = 2*q1;
    Hk2(2,4) = -2*q0;
    Hk2(3,1) = 2*q2;
    Hk2(3,2) = 2*q3;
    Hk2(3,3) = 2*q0;
    Hk2(3,4) = 2*q1;
    SK_MagX = zeros(1,1);
    SK_MagX(1) = 1/(R_Mag + 2*q0*(2*Pk_k1(1,1)*q0 + 2*Pk_k1(2,1)*q1 - 2*Pk_k1(3,1)*q2 - 2*Pk_k1(4,1)*q3) + 2*q1*(2*Pk_k1(1,2)*q0 + 2*Pk_k1(2,2)*q1 - 2*Pk_k1(3,2)*q2 - 2*Pk_k1(4,2)*q3) - 2*q2*(2*Pk_k1(1,3)*q0 + 2*Pk_k1(2,3)*q1 - 2*Pk_k1(3,3)*q2 - 2*Pk_k1(4,3)*q3) - 2*q3*(2*Pk_k1(1,4)*q0 + 2*Pk_k1(2,4)*q1 - 2*Pk_k1(3,4)*q2 - 2*Pk_k1(4,4)*q3));
    SK_MagY = zeros(1,1);
    SK_MagY(1) = 1/(R_Mag + 2*q3*(2*Pk_k1(1,1)*q3 - 2*Pk_k1(2,1)*q2 - 2*Pk_k1(3,1)*q1 + 2*Pk_k1(4,1)*q0) - 2*q2*(2*Pk_k1(1,2)*q3 - 2*Pk_k1(2,2)*q2 - 2*Pk_k1(3,2)*q1 + 2*Pk_k1(4,2)*q0) - 2*q1*(2*Pk_k1(1,3)*q3 - 2*Pk_k1(2,3)*q2 - 2*Pk_k1(3,3)*q1 + 2*Pk_k1(4,3)*q0) + 2*q0*(2*Pk_k1(1,4)*q3 - 2*Pk_k1(2,4)*q2 - 2*Pk_k1(3,4)*q1 + 2*Pk_k1(4,4)*q0));
    SK_MagZ = zeros(1,1);
    SK_MagZ(1) = 1/(R_Mag + 2*q2*(2*Pk_k1(1,1)*q2 + 2*Pk_k1(2,1)*q3 + 2*Pk_k1(3,1)*q0 + 2*Pk_k1(4,1)*q1) + 2*q3*(2*Pk_k1(1,2)*q2 + 2*Pk_k1(2,2)*q3 + 2*Pk_k1(3,2)*q0 + 2*Pk_k1(4,2)*q1) + 2*q0*(2*Pk_k1(1,3)*q2 + 2*Pk_k1(2,3)*q3 + 2*Pk_k1(3,3)*q0 + 2*Pk_k1(4,3)*q1) + 2*q1*(2*Pk_k1(1,4)*q2 + 2*Pk_k1(2,4)*q3 + 2*Pk_k1(3,4)*q0 + 2*Pk_k1(4,4)*q1));
    Kk2 = zeros(4,3);
    Kk2(1,1) = SK_MagX(1)*(2*Pk_k1(1,1)*q0 + 2*Pk_k1(1,2)*q1 - 2*Pk_k1(1,3)*q2 - 2*Pk_k1(1,4)*q3);
    Kk2(1,2) = -SK_MagY(1)*(2*Pk_k1(1,1)*q3 - 2*Pk_k1(1,2)*q2 - 2*Pk_k1(1,3)*q1 + 2*Pk_k1(1,4)*q0);
    Kk2(1,3) = SK_MagZ(1)*(2*Pk_k1(1,1)*q2 + 2*Pk_k1(1,3)*q0 + 2*Pk_k1(1,2)*q3 + 2*Pk_k1(1,4)*q1);
    Kk2(2,1) = SK_MagX(1)*(2*Pk_k1(2,1)*q0 + 2*Pk_k1(2,2)*q1 - 2*Pk_k1(2,3)*q2 - 2*Pk_k1(2,4)*q3);
    Kk2(2,2) = -SK_MagY(1)*(2*Pk_k1(2,1)*q3 - 2*Pk_k1(2,2)*q2 - 2*Pk_k1(2,3)*q1 + 2*Pk_k1(2,4)*q0);
    Kk2(2,3) = SK_MagZ(1)*(2*Pk_k1(2,1)*q2 + 2*Pk_k1(2,3)*q0 + 2*Pk_k1(2,2)*q3 + 2*Pk_k1(2,4)*q1);
    Kk2(3,1) = SK_MagX(1)*(2*Pk_k1(3,1)*q0 + 2*Pk_k1(3,2)*q1 - 2*Pk_k1(3,3)*q2 - 2*Pk_k1(3,4)*q3);
    Kk2(3,2) = -SK_MagY(1)*(2*Pk_k1(3,1)*q3 - 2*Pk_k1(3,2)*q2 - 2*Pk_k1(3,3)*q1 + 2*Pk_k1(3,4)*q0);
    Kk2(3,3) = SK_MagZ(1)*(2*Pk_k1(3,1)*q2 + 2*Pk_k1(3,3)*q0 + 2*Pk_k1(3,2)*q3 + 2*Pk_k1(3,4)*q1);
    Kk2(4,1) = SK_MagX(1)*(2*Pk_k1(4,1)*q0 + 2*Pk_k1(4,2)*q1 - 2*Pk_k1(4,3)*q2 - 2*Pk_k1(4,4)*q3);
    Kk2(4,2) = -SK_MagY(1)*(2*Pk_k1(4,1)*q3 - 2*Pk_k1(4,2)*q2 - 2*Pk_k1(4,3)*q1 + 2*Pk_k1(4,4)*q0);
    Kk2(4,3) = SK_MagZ(1)*(2*Pk_k1(4,1)*q2 + 2*Pk_k1(4,3)*q0 + 2*Pk_k1(4,2)*q3 + 2*Pk_k1(4,4)*q1);
    hk2 = zeros(3,1);
    hk2(1) = q0^2 + q1^2 - q2^2 - q3^2;
    hk2(2) = 2*q1*q2 - 2*q0*q3;
    hk2(3) = 2*q0*q2 + 2*q1*q3;
    Sqc2 = zeros(3,1);
    Sqc2(1) = hk2(3,1) - magZ;
    Sqc2(2) = hk2(2,1) - magY;
    Sqc2(3) = hk2(1,1) - magX;
    qc2 = zeros(4,1);
    qc2(1) = - Kk2(1,1)*Sqc2(3) - Kk2(1,2)*Sqc2(2) - Kk2(1,3)*Sqc2(1);
    qc2(2) = - Kk2(2,1)*Sqc2(3) - Kk2(2,2)*Sqc2(2) - Kk2(2,3)*Sqc2(1);
    qc2(3) = - Kk2(3,1)*Sqc2(3) - Kk2(3,2)*Sqc2(2) - Kk2(3,3)*Sqc2(1);
    qc2(4) = - Kk2(4,1)*Sqc2(3) - Kk2(4,2)*Sqc2(2) - Kk2(4,3)*Sqc2(1);
    SPk2 = zeros(16,1);
    SPk2(1) = Hk2(1,4)*Kk2(4,1) + Hk2(2,4)*Kk2(4,2) + Hk2(3,4)*Kk2(4,3);
    SPk2(2) = Hk2(1,3)*Kk2(4,1) + Hk2(2,3)*Kk2(4,2) + Hk2(3,3)*Kk2(4,3);
    SPk2(3) = Hk2(1,2)*Kk2(4,1) + Hk2(2,2)*Kk2(4,2) + Hk2(3,2)*Kk2(4,3);
    SPk2(4) = Hk2(1,1)*Kk2(4,1) + Hk2(2,1)*Kk2(4,2) + Hk2(3,1)*Kk2(4,3);
    SPk2(5) = Hk2(1,4)*Kk2(3,1) + Hk2(2,4)*Kk2(3,2) + Hk2(3,4)*Kk2(3,3);
    SPk2(6) = Hk2(1,3)*Kk2(3,1) + Hk2(2,3)*Kk2(3,2) + Hk2(3,3)*Kk2(3,3);
    SPk2(7) = Hk2(1,2)*Kk2(3,1) + Hk2(2,2)*Kk2(3,2) + Hk2(3,2)*Kk2(3,3);
    SPk2(8) = Hk2(1,1)*Kk2(3,1) + Hk2(2,1)*Kk2(3,2) + Hk2(3,1)*Kk2(3,3);
    SPk2(9) = Hk2(1,4)*Kk2(2,1) + Hk2(2,4)*Kk2(2,2) + Hk2(3,4)*Kk2(2,3);
    SPk2(10) = Hk2(1,3)*Kk2(2,1) + Hk2(2,3)*Kk2(2,2) + Hk2(3,3)*Kk2(2,3);
    SPk2(11) = Hk2(1,2)*Kk2(2,1) + Hk2(2,2)*Kk2(2,2) + Hk2(3,2)*Kk2(2,3);
    SPk2(12) = Hk2(1,1)*Kk2(2,1) + Hk2(2,1)*Kk2(2,2) + Hk2(3,1)*Kk2(2,3);
    SPk2(13) = Hk2(1,4)*Kk2(1,1) + Hk2(2,4)*Kk2(1,2) + Hk2(3,4)*Kk2(1,3);
    SPk2(14) = Hk2(1,3)*Kk2(1,1) + Hk2(2,3)*Kk2(1,2) + Hk2(3,3)*Kk2(1,3);
    SPk2(15) = Hk2(1,2)*Kk2(1,1) + Hk2(2,2)*Kk2(1,2) + Hk2(3,2)*Kk2(1,3);
    SPk2(16) = Hk2(1,1)*Kk2(1,1) + Hk2(2,1)*Kk2(1,2) + Hk2(3,1)*Kk2(1,3);
    Pk = zeros(4,4);
    Pk(1,1) = Pk1(1,1) - Pk1(1,1)*SPk2(16) - Pk1(2,1)*SPk2(15) - Pk1(3,1)*SPk2(14) - Pk1(4,1)*SPk2(13);
    Pk(1,2) = Pk1(1,2) - Pk1(1,2)*SPk2(16) - Pk1(2,2)*SPk2(15) - Pk1(3,2)*SPk2(14) - Pk1(4,2)*SPk2(13);
    Pk(1,3) = Pk1(1,3) - Pk1(1,3)*SPk2(16) - Pk1(2,3)*SPk2(15) - Pk1(3,3)*SPk2(14) - Pk1(4,3)*SPk2(13);
    Pk(1,4) = Pk1(1,4) - Pk1(1,4)*SPk2(16) - Pk1(2,4)*SPk2(15) - Pk1(3,4)*SPk2(14) - Pk1(4,4)*SPk2(13);
    Pk(2,1) = Pk1(2,1) - Pk1(1,1)*SPk2(12) - Pk1(2,1)*SPk2(11) - Pk1(3,1)*SPk2(10) - Pk1(4,1)*SPk2(9);
    Pk(2,2) = Pk1(2,2) - Pk1(1,2)*SPk2(12) - Pk1(2,2)*SPk2(11) - Pk1(3,2)*SPk2(10) - Pk1(4,2)*SPk2(9);
    Pk(2,3) = Pk1(2,3) - Pk1(1,3)*SPk2(12) - Pk1(2,3)*SPk2(11) - Pk1(3,3)*SPk2(10) - Pk1(4,3)*SPk2(9);
    Pk(2,4) = Pk1(2,4) - Pk1(1,4)*SPk2(12) - Pk1(2,4)*SPk2(11) - Pk1(3,4)*SPk2(10) - Pk1(4,4)*SPk2(9);
    Pk(3,1) = Pk1(3,1) - Pk1(1,1)*SPk2(8) - Pk1(2,1)*SPk2(7) - Pk1(3,1)*SPk2(6) - Pk1(4,1)*SPk2(5);
    Pk(3,2) = Pk1(3,2) - Pk1(1,2)*SPk2(8) - Pk1(2,2)*SPk2(7) - Pk1(3,2)*SPk2(6) - Pk1(4,2)*SPk2(5);
    Pk(3,3) = Pk1(3,3) - Pk1(1,3)*SPk2(8) - Pk1(2,3)*SPk2(7) - Pk1(3,3)*SPk2(6) - Pk1(4,3)*SPk2(5);
    Pk(3,4) = Pk1(3,4) - Pk1(1,4)*SPk2(8) - Pk1(2,4)*SPk2(7) - Pk1(3,4)*SPk2(6) - Pk1(4,4)*SPk2(5);
    Pk(4,1) = Pk1(4,1) - Pk1(1,1)*SPk2(4) - Pk1(2,1)*SPk2(3) - Pk1(3,1)*SPk2(2) - Pk1(4,1)*SPk2(1);
    Pk(4,2) = Pk1(4,2) - Pk1(1,2)*SPk2(4) - Pk1(2,2)*SPk2(3) - Pk1(3,2)*SPk2(2) - Pk1(4,2)*SPk2(1);
    Pk(4,3) = Pk1(4,3) - Pk1(1,3)*SPk2(4) - Pk1(2,3)*SPk2(3) - Pk1(3,3)*SPk2(2) - Pk1(4,3)*SPk2(1);
    Pk(4,4) = Pk1(4,4) - Pk1(1,4)*SPk2(4) - Pk1(2,4)*SPk2(3) - Pk1(3,4)*SPk2(2) - Pk1(4,4)*SPk2(1);
    
    qc1(4)=0;qc2(2)=0;qc2(3)=0;
    qk1=q+qc1;
    q=qk1+qc2;
    
    time=time+toc;   
end
time
a=(1:1:N);
figure(1);
subplot(3,1,1);
plot(a,Phi);
xlabel('k');ylabel('Phi');
subplot(3,1,2);
plot(a,IMU(:,10));
xlabel('k');ylabel('Phi_d');
subplot(3,1,3);
Phi_err=Phi'-IMU(:,10);
plot(a,Phi_err);
xlabel('k');ylabel('Phi_i');

figure(2);
subplot(3,1,1);
plot(a,Theta);
xlabel('k');ylabel('Theta');
subplot(3,1,2);
plot(a,IMU(:,11));
xlabel('k');ylabel('Theta_d');
subplot(3,1,3);
Theta_err=Theta'-IMU(:,11);
plot(a,Theta_err);
xlabel('k');ylabel('Theta_i');

figure(3);
subplot(3,1,1);
plot(a,Psi);
xlabel('k');ylabel('Psi');
subplot(3,1,2);
plot(a,IMU(:,12));
xlabel('k');ylabel('Psi_d');
subplot(3,1,3);
Psi_err=Psi'-IMU(:,12);
plot(a,Psi_err);
xlabel('k');ylabel('Psi_i');
