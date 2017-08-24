% IMPORTANT - This script requires the Matlab symbolic toolbox
% The code is the implementations of Double-Stage Kalman.

% Author: Xinzhe Gui, Jing Xue
% Email: lwft@qq.com
%%
clear;
close all;

addpath(genpath(pwd));

%加载数据
IMU=xlsread('TestData/datal_1.xlsx','Sheet2');

%my 6 points
AccSensx=0.995544148719357;
AccSensy=1.00067100989279;
AccSensz=0.999918078434263;
AccOffsetx=0.00428358112745599;
AccOffsety=-0.00270633319254957;
AccOffsetz=-0.0705557225725813;

magcalibration=[2.06385073217840,0.0413407668037351,0.0246229991758043;0,2.10288949976811,-0.0104899029262219;0,0,2.36140472565193];
magoffset=[0.0610603654152471;0.0201527841522371;0.990707291320690];

Gyro_x_offset=-19.602; Gyro_y_offset=11.886; Gyro_z_offset=1.338;

%初始值设定
g=9.8;unit=32768;

[mIMUData,nIMUData]=size(IMU);
N=mIMUData;
deltet=0.01;
Q=10^(-6)*eye(4,4);
Vk1=1;
% Rk1=18*8*g/unit*eye(3,3);
Rk1=2*eye(3,3);
Vk2=1;
% Rk2=1.0e-03*0.2327*eye(3,3);
Rk2=eye(3,3);

Mag_unit=0.3;
pi=3.141592653;
% q=angle2quat(IMU(1,12)*pi/180,IMU(1,11)*pi/180,IMU(1,10)*pi/180)';
% theta=90-59.0667;
theta=0;

q=[1;0;0;0];
q_1=q;
% Pk=eye(4,4)*(0.001-0.0002)+ones(4,4)*0.0002;
Pk=eye(4,4)*(0.125-0.0003)+ones(4,4)*0.0003;
Theta=zeros(1,N);
Phi=zeros(1,N);
Psi=zeros(1,N);
Theta2=zeros(1,N);
Phi2=zeros(1,N);
Psi2=zeros(1,N);
Cbn_save(3,3,N)=0;
time=0;
for k=1:N

  normq=sqrt(q(1)^2+q(2)^2+q(3)^2+q(4)^2);
  q(1)=q(1)/normq;
  q(2)=q(2)/normq;
  q(3)=q(3)/normq;
  q(4)=q(4)/normq;
  
  normq_1=sqrt(q_1(1)^2+q_1(2)^2+q_1(3)^2+q_1(4)^2);
  q_1(1)=q_1(1)/normq_1;
  q_1(2)=q_1(2)/normq_1;
  q_1(3)=q_1(3)/normq_1;
  q_1(4)=q_1(4)/normq_1;
   %方向旋转矩阵
  Cnb = [q(1)^2+q(2)^2-(q(3)^2+q(4)^2),     2*(q(2)*q(3)+q(1)*q(4)),     2*(q(2)*q(4)-q(1)*q(3));
       2*(q(2)*q(3)-q(1)*q(4)),         q(1)^2+q(3)^2-(q(2)^2+q(4)^2),     2*(q(1)*q(2)+q(3)*q(4));
       2*(q(2)*q(4)+q(1)*q(3)),     2*(q(3)*q(4)-q(1)*q(2)),        q(1)^2+q(4)^2-(q(2)^2+q(3)^2)];

   Cbn=Cnb';
   Cbn_save(:,:,k)=Cbn;
  
  Phi(k)=atan2(2*(q(3)*q(4)+q(1)*q(2)),q(1)^2+q(4)^2-(q(2)^2+q(3)^2))/pi*180;%ROLL*pi/180  
  Theta(k)=asin(2*(q(1)*q(3)-q(4)*q(2)))/pi*180;%PITCH
  Psi(k)= atan2(2*(q(2)*q(3)+q(1)*q(4)),(q(1)^2+q(2)^2-q(3)^2-q(4)^2))/pi*180;%YAW

  Phi2(k)=atan2(2*(q_1(3)*q_1(4)+q_1(1)*q_1(2)),q_1(1)^2+q_1(4)^2-(q_1(2)^2+q_1(3)^2))/pi*180;%ROLL*pi/180  
  Theta2(k)=asin(2*(q_1(1)*q_1(3)-q_1(4)*q_1(2)))/pi*180;%PITCH
  Psi2(k)= atan2(2*(q_1(2)*q_1(3)+q_1(1)*q_1(4)),(q_1(1)^2+q_1(2)^2-q_1(3)^2-q_1(4)^2))/pi*180;%YAW

  tic
  %读取陀螺仪的值
  wx=(IMU(k,4)-Gyro_x_offset)*2000/unit*pi/180;
  wy=-(IMU(k,5)-Gyro_y_offset)*2000/unit*pi/180;
  wz=-(IMU(k,6)-Gyro_z_offset)*2000/unit*pi/180;
  deltet=IMU(k,17)/1000000;
  
  %旋转矩阵的变化速率
  Omega=[0  -wx -wy -wz;
         wx  0   wz -wy;
         wy -wz  0   wx;
         wz  wy -wx  0 ];
  A=eye(4,4)+1/2*Omega*deltet;
  q=A*q;
  q_1=A*q_1;

  Pk_k1=A*Pk*A'+Q;
  %stage 1
%   Hk1=[-2*q(3)  2*q(4) -2*q(1) 2*q(2);
%         2*q(2)  2*q(1)  2*q(4) 2*q(3);
%         2*q(1) -2*q(2) -2*q(3) 2*q(4)];
  Hk1=g*[-2*q(3)  2*q(4) -2*q(1) 2*q(2);
        2*q(2)  2*q(1)  2*q(4) 2*q(3);
        2*q(1) -2*q(2) -2*q(3) 2*q(4)];    
   Kk1=Pk_k1*Hk1'*(Hk1*Pk_k1*Hk1'+Vk1*Rk1*Vk1')^(-1);
   Kk1_X=Pk_k1*Hk1(1,:)'*(Hk1(1,:)*Pk_k1*Hk1(1,:)'+2)^(-1);
   Kk1_Y=Pk_k1*Hk1(2,:)'*(Hk1(2,:)*Pk_k1*Hk1(2,:)'+2)^(-1);
   Kk1_Z=Pk_k1*Hk1(3,:)'*(Hk1(3,:)*Pk_k1*Hk1(3,:)'+2)^(-1);
   Kk1_temp=[Kk1_X Kk1_Y Kk1_Z];
   Kk1=Kk1_temp;
   %读取加速度的值
   amx=-((IMU(k,1)*8/unit)-AccOffsetx)*AccSensx*g;
   amy=((IMU(k,2)*8/unit)-AccOffsety)*AccSensy*g;
   amz=((IMU(k,3)*8/unit)-AccOffsetz)*AccSensz*g;
   
   UKU(k,:)=Cbn*[amx;amy;amz];
  
   Zk1=[amx;amy;amz];
   
   hk1=g*[2*(q(2)*q(4)-q(1)*q(3));2*(q(1)*q(2)+q(3)*q(4));q(1)^2+q(4)^2-(q(2)^2+q(3)^2)];

   qc1=Kk1*(Zk1-hk1);
   qc1(4)=0;
   qk1=q+qc1;
   Pk1=(eye(4,4)-Kk1*Hk1)*Pk_k1;
   
   q1=q(1);
   q2=q(2);
   q3=q(3);
   q4=q(4);
   
   %stage 2
   Hk2 = [ 2*q1*cos(theta) - 2*q3*sin(theta), 2*q2*cos(theta) + 2*q4*sin(theta), - 2*q3*cos(theta) - 2*q1*sin(theta), 2*q2*sin(theta) - 2*q4*cos(theta);
       2*q2*sin(theta) - 2*q4*cos(theta), 2*q3*cos(theta) + 2*q1*sin(theta),   2*q2*cos(theta) + 2*q4*sin(theta), 2*q3*sin(theta) - 2*q1*cos(theta);
       2*q3*cos(theta) + 2*q1*sin(theta), 2*q4*cos(theta) - 2*q2*sin(theta),   2*q1*cos(theta) - 2*q3*sin(theta), 2*q2*cos(theta) + 2*q4*sin(theta)];

   Kk2=Pk_k1*Hk2'*(Hk2*Pk_k1*Hk2'+Vk2*Rk2*Vk2')^(-1);
   Kk2_X=Pk_k1*Hk2(1,:)'*(Hk2(1,:)*Pk_k1*Hk2(1,:)'+1)^(-1);
   Kk2_Y=Pk_k1*Hk2(2,:)'*(Hk2(2,:)*Pk_k1*Hk2(2,:)'+1)^(-1);
   Kk2_Z=Pk_k1*Hk2(3,:)'*(Hk2(3,:)*Pk_k1*Hk2(3,:)'+1)^(-1);
   Kk2_temp=[Kk2_X Kk2_Y Kk2_Z];
   Kk2=Kk2_temp;
   
   %读取磁力计的值
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

   Zk2=[mx_s;my_s;mz_s];

   hk2=[cos(theta)*(q1^2 + q2^2 - q3^2 - q4^2) - sin(theta)*(2*q1*q3 - 2*q2*q4);
       sin(theta)*(2*q1*q2 + 2*q3*q4) - cos(theta)*(2*q1*q4 - 2*q2*q3);
       cos(theta)*(2*q1*q3 + 2*q2*q4) + sin(theta)*(q1^2 - q2^2 - q3^2 + q4^2)];

   qc2=Kk2*(Zk2-hk2);

   qc2(2)=0;qc2(3)=0;
   q=qk1+qc2;
   Pk=(eye(4,4)-Kk2*Hk2)*Pk1;
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
xlabel('k');ylabel('Phi_err');

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
xlabel('k');ylabel('Theta_err');

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
xlabel('k');ylabel('Psi_err');
     