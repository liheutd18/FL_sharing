function  [YPred,YTest,net,sig,mu,MAPE]=LSTM_1(TotalL,L_dis,Newl,Newntl,TotalNTL,Totalestimatepv,estimatepv,stdpv)%N*24
    for j=1:5
%35040 points. select start and end, 1day=96 points
startday=1;%75 168 280
d=30;%104 197 309
data = TotalNTL((startday-1)*24+1:24*d)';
%data = Newntl((startday-1)*24+1:24*d,1)'; %Newntl only or L_dis
data_stdpv=stdpv((startday-1)*24+1:24*d);
data=[data ; data_stdpv];
   %data_pv = Totalestimatepv((startday-1)*24+1:24*d,1)'; % estimate pv or Totalestimatepv
%data_realwobtm = Newl((startday-1)*24+1:24*d,1)';
%    data_realntl = Newntl((startday-1)*24+1:24*d,1:2)'; %Newntl or TotalNTL
%    data_realntl=[data_realntl;Newntl((startday-1)*24+1:24*d,6:7)'];% ;Newntl((startday-1)*24+1:24*d,9)'];
%    data=[data ; data_realntl];


% data_corelation=[data ;Newntl((startday-1)*24+1:24*d,:)'];
% corrplot(data_corelation(:,1:end-72)')
[row,col] =size(data);
dayindex=[1:24];K1=repmat(dayindex,1,365);
int=fix(d/7);rem=d-7*int;
%input int first day is monday in 2018
weeknum=[1:7];week7index=repmat(weeknum,1,54);
K2=[week7index, 1:rem];K2=ones(24,1)*K2;K2=K2(:);K2=K2';K2=K2(1:8760);

%holi=holidays('jan 1 2018', 'dec 31 2018');
k_holiday=[1 15 45 70 91 92 133 148 168 185 246 254 304 326 327 328 358 359 360 361 362 363 364 365];
K3=zeros(1,365);
    for i=1:length(k_holiday)
        K3(k_holiday(i))=1;
    end
K3=ones(24,1)*K3;K3=K3(:);K3=K3';

K=[K1;K2;K3];
K=K(:,(startday-1)*24+1:24*d);
figure
plot(data')
xlabel("time (1h)")
ylabel("load")
title("Netload forecasting")

        numTimeStepsTrain = floor(0.9*numel(data(1,:)))-1;
%        dataTrain = data(1:numTimeStepsTrain+1);

%        data=[data;K];
        dataTrain = data(:,1:numTimeStepsTrain+1); %TRAIN
        
        dataTest = data(:,numTimeStepsTrain+1:end); %TEST
%        dataTest_realwobtm = data_realwobtm(:,numTimeStepsTrain+1:end); %TEST
        
        numFeatures = row; %
        numResponses = row; %        
        mu=zeros(numFeatures,1);sig=zeros(numFeatures,1);
        for i=1:numFeatures
            mu(i) = mean(dataTrain(i,:));
            sig(i) = std(dataTrain(i,:));
        end
        dataTrainStandardized = (dataTrain - mu) ./ sig;
        XTrain = dataTrainStandardized(:,1:end-1);%remain1
        YTrain = dataTrainStandardized(:,2:end);% emit first one

        numHiddenUnits = 200;

        layers = [ ...
            sequenceInputLayer(numFeatures)
            lstmLayer(numHiddenUnits)
            fullyConnectedLayer(numResponses)
            regressionLayer];
        options = trainingOptions('adam', ...
            'MaxEpochs',150, ...
            'GradientThreshold',1, ...
            'InitialLearnRate',0.005, ...
            'LearnRateSchedule','piecewise', ...
            'LearnRateDropPeriod',125, ...
            'LearnRateDropFactor',0.2, ...
            'Verbose',0, ...
            'Plots','training-progress');
        net(j) = trainNetwork(XTrain,YTrain,layers,options);
%        load lstmnet0213
        dataTestStandardized = (dataTest - mu) ./ sig;
        XTest = dataTestStandardized(:,1:end-1);
        net(j) = predictAndUpdateState(net(j),XTrain);%update last point
        [net(j),YPred] = predictAndUpdateState(net(j),YTrain(:,end));

        numTimeStepsTest = numel(XTest(1,:));
        for i = 2:numTimeStepsTest
            [net(j),YPred(:,i)] = predictAndUpdateState(net(j),YPred(:,i-1),'ExecutionEnvironment','cpu');
        end
        YPred = sig.*YPred + mu;
        YPred_DA=YPred;
        YTest = dataTest(:,2:end);
%        YTest1 = dataTest_realwobtm(:,2:end);

        rmse1 = sqrt(mean((YPred(1,:)-YTest(1,:)).^2));
        nrmse1=rmse1/(max(YTest(1,:))-min(YTest(1,:)));
        
%         rmse11 = sqrt(mean((YPred(1,:)+data_pv(end-70:end)-data_realntl(end-70:end)).^2));
%         nrmse11=rmse11/(max(data_realntl(end-70:end))-min(data_realntl(end-70:end)));
        
        figure
        plot(dataTrain(1,1:end-1))
        hold on
        idx = numTimeStepsTrain:(numTimeStepsTrain+numTimeStepsTest);
        plot(idx,[data(1,numTimeStepsTrain) YPred(1,:)],'.-')
        hold off
        xlabel("time")
        ylabel("power kw")
        title("Forecast")
        legend(["Observed" "Forecast"])
        figure
        subplot(2,1,1)
        plot(YTest(1,:))
        hold on
        plot(YPred(1,:),'.-')
        hold off
        legend(["Observed" "Forecast"])
        ylabel("Cases")
        title("Forecast")
        hold on

        %resetnet
        net(j) = resetState(net(j));
        net(j) = predictAndUpdateState(net(j),XTrain);
        YPred = [];
        numTimeStepsTest = numel(XTest(1,:));
        XTest(1,:)=XTest(1,:).*(0.8+0.2*rand(numTimeStepsTest,1))';
        for i = 1:numTimeStepsTest
            [net(j),YPred(:,i)] = predictAndUpdateState(net(j),XTest(:,i),'ExecutionEnvironment','cpu');
        end
        YPred = sig.*YPred + mu;
        rmse2 = sqrt(mean((YPred(1,:)-YTest(1,:)).^2));
        nrmse2=rmse2/(max(YTest(1,:))-min(YTest(1,:)));
        
%         
%         rmse3 = sqrt(mean((YPred(1,:)-YTest1(1,:)).^2));
%         nrmse3=rmse2/(max(YTest1(1,:))-min(YTest1(1,:)));
        
%         rmse22 = sqrt(mean((YPred(1,:)+data_pv(end-70:end)-data_realntl(end-70:end)).^2));
%         nrmse22=rmse22/(max(data_realntl(end-70:end))-min(data_realntl(end-70:end)));
%%        

        subplot(2,1,2)
        plot(YTest(1,:))
        hold on
        plot(YPred(1,:),'.-')
        hold off
        legend(["Observed" "Predicted"])
        ylabel("Cases")
        title("Forecast with Updates")
        Results{j}={net(j),YPred,YPred_DA,YTest,nrmse1,nrmse2};   
        RMSE_all(j,:)=[nrmse1 nrmse2];
    end
end


