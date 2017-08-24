% IMPORTANT - This script requires the Matlab symbolic toolbox
% The code is the implementations of Double-Stage Kalman.

% Author: Xinzhe Gui
% Email: lwft@qq.com
%% environment config
clear all;close all;
reset(symengine);
addpath(genpath(pwd));
%%
% this code generation use the CovariancePrediction_EasyQ.m file, not the
% CovariancePrediction.
% the syms behind is the necessary syms for double-stage kalman.
% the symbol "//" before the syms means that they are used for
% CovariancePrediction, not the CovariancePrediction_EasyQ.

% quatNew
% SF
% //SQ
% //Q
% //SG
% //G
% Pk_k1
% 
% Hk1
% SK_dVelX
% SK_dVelY
% SK_dVelZ
% Kk1
% hk1
% Sqc1
% qc1
% SPk1
% Pk1
% 
% Hk2
% SK_MagX
% SK_MagY
% SK_MagZ
% Kk2
% hk2
% Sqc2
% qc2
% SPk2
% Pk2
%% Load Data
load('./Mat/StatePrediction.mat');
% load('./Mat/CovariancePrediction.mat');
load('./Mat/CovariancePrediction_EasyQ.mat');
load('./Mat/dVel.mat');
load('./Mat/Mag.mat');
load('./Mat/qc1AndPk1.mat');
load('./Mat/qc2AndPk2.mat');

%% Set the file address and Open output file
SourceCodeFileName = strcat('./GeneratedCode/','SourceCode.txt');
MatlabCodeFileName = strcat('./GeneratedCode/','MatlabCode.m');
CCodeFileName = strcat('./GeneratedCode/','CCode.c');

fileName = SourceCodeFileName;
fid = fopen(fileName,'wt');
%% Write equation for one step state prediction
if exist('quatNew','var')
    
    fprintf(fid,'quatNew = zeros(%d,1);\n',numel(quatNew));
    for rowIndex = 1:numel(quatNew)
        string = char(quatNew(rowIndex,1));
        fprintf(fid,'quatNew(%d) = %s;\n',rowIndex,string);
    end  
    fprintf(fid,'\n');
    
end
%% Write equation for state transition matrix
if exist('SF','var')
    
    fprintf(fid,'SF = zeros(%d,1);\n',numel(SF));
    for rowIndex = 1:numel(SF)
        string = char(SF(rowIndex,1));
        fprintf(fid,'SF(%d) = %s;\n',rowIndex,string);
    end  
    fprintf(fid,'\n');
    
end
%% Write equations for covariance prediction
% if exist('SPk_k1','var')
%     
%     fprintf(fid,'SPk_k1 = zeros(%d,1);\n',numel(SPk_k1));
%     for rowIndex = 1:numel(SPk_k1)
%         string = char(SPk_k1(rowIndex,1));
%         fprintf(fid,'SPk_k1(%d) = %s;\n',rowIndex,string);
%     end
%     fprintf(fid,'\n');
%     
% end
%% Write equations for covariance prediction
if exist('Pk_k1','var')
    
    [nRow,nCol] = size(Pk_k1);
    fprintf(fid,'Pk_k1 = zeros(%d,%d);\n',nRow,nCol);
    for rowIndex = 1:nRow
        for colIndex = 1:nCol
            string = char(Pk_k1(rowIndex,colIndex));
            % don't write out a zero-assignment
            if ~strcmpi(string,'0')
                fprintf(fid,'Pk_k1(%d,%d) = %s;\n',rowIndex,colIndex,string);
            end
        end
    end
    fprintf(fid,'\n');
    
end
%% Write equations for accelerator data fusion
if exist('Hk1','var')
    
    [nRow,nCol] = size(Hk1);
    fprintf(fid,'Hk1 = zeros(%d,%d);\n',nRow,nCol);
    for rowIndex = 1:nRow
        for colIndex = 1:nCol
            string = char(Hk1(rowIndex,colIndex));
            % don't write out a zero-assignment
            if ~strcmpi(string,'0')
                fprintf(fid,'Hk1(%d,%d) = %s;\n',rowIndex,colIndex,string);
            end
        end
    end
    fprintf(fid,'\n');
    
end
%% Write Kalman gain equations for accelerator data fusion
if exist('SK_dVelX','var')
    
    fprintf(fid,'SK_dVelX = zeros(%d,1);\n',numel(SK_dVelX));
    for rowIndex = 1:numel(SK_dVelX)
        string = char(SK_dVelX(rowIndex,1));
        fprintf(fid,'SK_dVelX(%d) = %s;\n',rowIndex,string);
    end
    fprintf(fid,'\n');
    
end
if exist('SK_dVelY','var')
    
    fprintf(fid,'SK_dVelY = zeros(%d,1);\n',numel(SK_dVelY));
    for rowIndex = 1:numel(SK_dVelY)
        string = char(SK_dVelY(rowIndex,1));
        fprintf(fid,'SK_dVelY(%d) = %s;\n',rowIndex,string);
    end
    fprintf(fid,'\n');
    
end
if exist('SK_dVelZ','var')
    
    fprintf(fid,'SK_dVelZ = zeros(%d,1);\n',numel(SK_dVelZ));
    for rowIndex = 1:numel(SK_dVelZ)
        string = char(SK_dVelZ(rowIndex,1));
        fprintf(fid,'SK_dVelZ(%d) = %s;\n',rowIndex,string);
    end
    fprintf(fid,'\n'); 
    
end
Kk1=[K_dVelX K_dVelY K_dVelZ];

if exist('Kk1','var')
    
    [nRow,nCol] = size(Kk1);
    fprintf(fid,'Kk1 = zeros(%d,%d);\n',nRow,nCol);
    for rowIndex = 1:nRow
        for colIndex = 1:nCol
            string = char(Kk1(rowIndex,colIndex));
            % don't write out a zero-assignment
            if ~strcmpi(string,'0')
                fprintf(fid,'Kk1(%d,%d) = %s;\n',rowIndex,colIndex,string);
            end
        end
    end
    fprintf(fid,'\n');
    
end
%% Write observation data for accelerator data fusion
if exist('hk1','var')
    
    fprintf(fid,'hk1 = zeros(%d,1);\n',numel(hk1));
    for rowIndex = 1:numel(hk1)
        string = char(hk1(rowIndex,1));
        fprintf(fid,'hk1(%d) = %s;\n',rowIndex,string);
    end
    fprintf(fid,'\n');
    
end
%% Write correction equations for accelerator data fusion
if exist('Sqc1','var')
    
    fprintf(fid,'Sqc1 = zeros(%d,1);\n',numel(Sqc1));
    for rowIndex = 1:numel(Sqc1)
        string = char(Sqc1(rowIndex,1));
        fprintf(fid,'Sqc1(%d) = %s;\n',rowIndex,string);
    end
    fprintf(fid,'\n');
    
end

if exist('qc1','var')
    
    fprintf(fid,'qc1 = zeros(%d,1);\n',numel(qc1));
    for rowIndex = 1:numel(qc1)
        string = char(qc1(rowIndex,1));
        fprintf(fid,'qc1(%d) = %s;\n',rowIndex,string);
    end
    fprintf(fid,'\n');
    
end
%% Write covariance equations for accelerator data fusion
if exist('SPk1','var')
    
    fprintf(fid,'SPk1 = zeros(%d,1);\n',numel(SPk1));
    for rowIndex = 1:numel(SPk1)
        string = char(SPk1(rowIndex,1));
        fprintf(fid,'SPk1(%d) = %s;\n',rowIndex,string);
    end
    fprintf(fid,'\n');
    
end

if exist('Pk1','var')
    
    [nRow,nCol] = size(Pk1);
    fprintf(fid,'Pk1 = zeros(%d,%d);\n',nRow,nCol);
    for rowIndex = 1:nRow
        for colIndex = 1:nCol
            string = char(Pk1(rowIndex,colIndex));
            % don't write out a zero-assignment
            if ~strcmpi(string,'0')
                fprintf(fid,'Pk1(%d,%d) = %s;\n',rowIndex,colIndex,string);
            end
        end
    end
    fprintf(fid,'\n');
    
end
%% Write equations for magnetometer data fusion
if exist('Hk2','var')
    
    [nRow,nCol] = size(Hk2);
    fprintf(fid,'Hk2 = zeros(%d,%d);\n',nRow,nCol);
    for rowIndex = 1:nRow
        for colIndex = 1:nCol
            string = char(Hk2(rowIndex,colIndex));
            % don't write out a zero-assignment
            if ~strcmpi(string,'0')
                fprintf(fid,'Hk2(%d,%d) = %s;\n',rowIndex,colIndex,string);
            end
        end
    end
    fprintf(fid,'\n');
    
end
%% Write Kalman gain equations for magnetometer data fusion
if exist('SK_MagX','var')
    
    fprintf(fid,'SK_MagX = zeros(%d,1);\n',numel(SK_MagX));
    for rowIndex = 1:numel(SK_MagX)
        string = char(SK_MagX(rowIndex,1));
        fprintf(fid,'SK_MagX(%d) = %s;\n',rowIndex,string);
    end
    fprintf(fid,'\n');
    
end
if exist('SK_MagY','var')
    
    fprintf(fid,'SK_MagY = zeros(%d,1);\n',numel(SK_MagY));
    for rowIndex = 1:numel(SK_MagY)
        string = char(SK_MagY(rowIndex,1));
        fprintf(fid,'SK_MagY(%d) = %s;\n',rowIndex,string);
    end
    fprintf(fid,'\n');
    
end
if exist('SK_MagZ','var')
    
    fprintf(fid,'SK_MagZ = zeros(%d,1);\n',numel(SK_MagZ));
    for rowIndex = 1:numel(SK_MagZ)
        string = char(SK_MagZ(rowIndex,1));
        fprintf(fid,'SK_MagZ(%d) = %s;\n',rowIndex,string);
    end
    fprintf(fid,'\n'); 
    
end
Kk2=[K_MagX K_MagY K_MagZ];

if exist('Kk2','var')
    
    [nRow,nCol] = size(Kk2);
    fprintf(fid,'Kk2 = zeros(%d,%d);\n',nRow,nCol);
    for rowIndex = 1:nRow
        for colIndex = 1:nCol
            string = char(Kk2(rowIndex,colIndex));
            % don't write out a zero-assignment
            if ~strcmpi(string,'0')
                fprintf(fid,'Kk2(%d,%d) = %s;\n',rowIndex,colIndex,string);
            end
        end
    end
    fprintf(fid,'\n');
    
end
%% Write observation data for magnetometer data fusion
if exist('hk2','var')
    
    fprintf(fid,'hk2 = zeros(%d,1);\n',numel(hk2));
    for rowIndex = 1:numel(hk2)
        string = char(hk2(rowIndex,1));
        fprintf(fid,'hk2(%d) = %s;\n',rowIndex,string);
    end
    fprintf(fid,'\n');
    
end
%% Write Correction equations for magnetometer data fusion
if exist('Sqc2','var')
    
    fprintf(fid,'Sqc2 = zeros(%d,1);\n',numel(Sqc2));
    for rowIndex = 1:numel(Sqc2)
        string = char(Sqc2(rowIndex,1));
        fprintf(fid,'Sqc2(%d) = %s;\n',rowIndex,string);
    end
    fprintf(fid,'\n');
    
end

if exist('qc2','var')
    
    fprintf(fid,'qc2 = zeros(%d,1);\n',numel(qc2));
    for rowIndex = 1:numel(qc2)
        string = char(qc2(rowIndex,1));
        fprintf(fid,'qc2(%d) = %s;\n',rowIndex,string);
    end
    fprintf(fid,'\n');
    
end
%% Write covariance equations for magnetometer data fusion
if exist('SPk2','var')
    
    fprintf(fid,'SPk2 = zeros(%d,1);\n',numel(SPk2));
    for rowIndex = 1:numel(SPk2)
        string = char(SPk2(rowIndex,1));
        fprintf(fid,'SPk2(%d) = %s;\n',rowIndex,string);
    end
    fprintf(fid,'\n');
    
end

if exist('Pk2','var')
    
    [nRow,nCol] = size(Pk2);
    fprintf(fid,'Pk2 = zeros(%d,%d);\n',nRow,nCol);
    for rowIndex = 1:nRow
        for colIndex = 1:nCol
            string = char(Pk2(rowIndex,colIndex));
            % don't write out a zero-assignment
            if ~strcmpi(string,'0')
                fprintf(fid,'Pk2(%d,%d) = %s;\n',rowIndex,colIndex,string);
            end
        end
    end
    fprintf(fid,'\n');
    
end
%% Close output file
fclose(fid);

ConvertToM(SourceCodeFileName,MatlabCodeFileName);
ConvertToC(SourceCodeFileName,CCodeFileName);