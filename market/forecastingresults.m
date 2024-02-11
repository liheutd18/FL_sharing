%load forecasting results
clc
clear
     Data=readtable('Forecasting results.xlsx','Sheet','sheet3');
     for j=1:10
         NetLoad_real(j,:)=table2array(Data(2+(j-1)*4,50:73)); 
         NetLoad_pred(j,:)=table2array(Data(1+(j-1)*4,50:73)); 
     end
         AGGNetLoad_predDA=table2array(Data(42,50:73)); 
         AGGNetLoad_real=table2array(Data(43,50:73));      
         AGGNetLoad_predHA=table2array(Data(46,50:73)); 
     a=1;