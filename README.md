## The implementations of Double-Stage Kalman ##

The code is the implementations of Double-Stage Kalman, and the theory is from the paper named "A Double-Stage Kalman Filter for Orientation Tracking With an Integrated Processor in 9-D IMU"

  * Support Website: http://www.lowfet.com

### Directory Explanation

  * GeneratedCode: the SourceCode, MatlabCode and CCode generated by Main_8GenerateCandMatlabCode.m
  * Mat: the data files generated by seven 'Main' files( Main_1StatePrediction.m,...,and Main_7Qc2andPk2.m)
  * TestData: the test data colleted by the author(Xinzhe Gui)
  * Tools: the necessary tool files, which are from PX4 and Ardupilot

### File Explanation

  * Main_1StatePrediction.m: generate the necessary file(StatePrediction.mat), which is the base of the others
  * Main_2CovariancePrediction.m: generate Pk_k1 with a new method, which contains control matrix.
  * Main_3CovariancePrediction_EasyQ.m: generate Pk_k1 only with Q, the method is similar as the paper mentioned.(A Double-Stage Kalman Filter for Orientation Tracking With an Integrated Processor in 9-D IMU)
  * Main_4AccKMatrix.m: generate Kk1
  * Main_5Qc1andPk1.m: generate qc1 and Pk1
  * Main_6MagKMatrix.m: generate Kk2
  * Main_7Qc2andPk2.m: generate qc2 and Pk2(it's also the Pk)
  * Main_8GenerateCandMatlabCode.m: generate SourceCode, MatlabCode and CCode
  * Test_OriginalDoubleStageKalman.m: the matlab implementations of Double-Stage Kalman
  * Test_OptimisedDoubleStageKalman.m: the optimised matlab implementations of Double-Stage Kalman
  * SourceCode.txt: the source syms
  * MatlabCode.m: the optimised matlab code
  * CCode.c: the optimised c code
  
### Using
  * First: Execute the file(Main_1StatePrediction.m,...,and Main_7Qc2andPk2.m) by order, then you can get the data file(.mat) from the Mat directory
  * Second: Execute the Main_8GenerateCandMatlabCode.m, and you can get the SourceCode, MatlabCode and CCode.
  * Third: If you want to see the result of the code, you should Execute the Test_OriginalDoubleStageKalman.m and Test_OptimisedDoubleStageKalman.m.
  * Fourth: If you want to use the optimised matlab code, you can get it from MatlabCode.m, and you can set qc1(4)=0,qc2(2)=0,qc2(3)=0 to simply the calculate.
  * Fifth: The CCode lack the code, which is used to do Quaternion normalizing, so if you use the CCode, you should appendix it.Meanwhile, you should change the Pk2 into Pk and you can set qc1[3]=0,qc2[1]=0,qc2[2]=0 to simply the calculate
  
### Contributors

  * Xinzhe Gui
  * Jing Xue
  
### Contact Information

  * Email: lwft@qq.com