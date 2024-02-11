%data processing
clear all 
clc
%Austin_data=readtable('C:\Users\lxh180005\Desktop\APPLYSI\15minute_data_austin\15minute_data_austin.csv');
load Pecanstreetdata.mat
Austin_data=[Austin_data(:,1) Austin_data(:,2) Austin_data(:,32) Austin_data(:,68)];
Austin_data.TIME = datetime(Austin_data.local_15min, 'InputFormat', 'yyyy-MM-dd HH:mm:ssXXX', 'TimeZone', 'UTC-6:00');

length(Austin_data.dataid(Austin_data.dataid==661));
B=unique(Austin_data.dataid);

for i=1:length(B)
    L=Austin_data.dataid(Austin_data.dataid==B(i));
    L1(:,i)=Austin_data.grid(Austin_data.dataid==B(i));
    P1(:,i)=Austin_data.solar(Austin_data.dataid==B(i));
    L3=Austin_data.TIME(Austin_data.dataid==B(i));
end
BB = sortrows(Austin_data,'TIME');
BBB=BB(:,3:5);


%     fname1 = ['MonthEnergy',num2str(Year),num2str(MaxPeakMM(Year))];
%     %load(fname1);
%     eval(['Power','=','MonthEnergy',num2str(MaxPeakMM),';']);