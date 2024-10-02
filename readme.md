### Source Codes for "A theoretical model for compressible bubble dynamics considering phase transition and migration"


#### 1. Purpose
This is a C/C++ code for the theoretical model for compressible bubble dynamics considering phase transition and migration. It is used to compute the bubble dynamics in different circumstances considering phase transition and migration. 


#### 2. Usage
The code is written with C/C++. Thus, a C/C++ compiler is required. We recommend installing C/C++ support in Windows Visual Studio. 

#### 3. Contents

##### 1) The input files for the programme consist of three documents: 

- "case.dat", enter fundamental scenario details such as the computation time; 
- "constants.dat", enter informatio related to the properties of the flow field, such as the gravity acceleration; 
- "thermodynamic.dat", enter the parameters related to phase transitions at the bubble surface, such as the thermal conductivity coefficient.

(The definitions of parameters could be seen at the end of the parameters in each file)


##### 2) The output files for the programme consist of two parts:

- "bubble.dat", document the time evolution of bubble radius and displacement (If N bubbles are computed, there are N output files: bubble1.dat, bubble2.dat, ... , bubbleN.dat.); 
- "pressure.dat", document the time evolution of the pressure at the flow field; 




#### 3. References
Users can refer to the following paper for detailed theory and equations:

[A-Man Zhang,  Shi-Min Li,  Run-Ze Xu, Shao-Cong Pei, Shuai Li, and  Yun-Long Liu, A theoretical model for compressible bubble dynamics considering phase transition and migration. arXiv preprint arXiv:2410.00335.]
(https://arxiv.org/abs/2410.00335)

###### Should you have any questions regarding the programme, please do not hesitate to get in touch at lishimien@126.com.



