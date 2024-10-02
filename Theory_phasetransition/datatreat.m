hold on
AA=load('bubble1.dat');
plot(AA(:,1)*1e3+0,AA(:,2)*1e3)
hold on
plot(AA(:,1)*1e3,AA(:,3)*1e3)

hold on
AA=load('bubble2.dat');
plot(AA(:,1)*1e0+0,AA(:,2)*1e0)
hold on
plot(AA(:,1),AA(:,3))

AA=load('pressure.dat');
plot(AA(:,1),AA(:,2))

hold on
AA=load('SDS1204X_HD_Binary_C1_6_Analog_Trace.csv');
plot(AA(:,1),AA(:,2))

AA=importdata('RU.txt');
plot(AA.data(:,1),AA.data(:,2))