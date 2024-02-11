%function  [YPred,YTest,net,sig,mu,MAPE]=LSTM_2(TotalNTL)%N*24
clc
clear all
for k=1:10
    load pecanYL.mat
    for i=1:20
        NTL(:,i)=PecanYL(:,2*i-1);
        PV(:,i)=PecanYL(:,2*i);
    end
    NTL=NTL(:,5:15);
    PV=PV(:,5:15);
    NTL(:,6)=[];PV(:,6)=[];
%35040 points. select start and end, 1day=96 points;
% dayindex=[1:24];K1=repmat(dayindex,1,365);
% int=fix(d/7);rem=d-7*int;
% %input int first day is monday in 2018
% weeknum=[1:7];week7index=repmat(weeknum,1,54);
% K2=[week7index, 1:rem];K2=ones(24,1)*K2;K2=K2(:);K2=K2';K2=K2(1:8760);
% 
% %holi=holidays('jan 1 2018', 'dec 31 2018');
% k_holiday=[1 15 45 70 91 92 133 148 168 185 246 254 304 326 327 328 358 359 360 361 362 363 364 365];
% K3=zeros(1,365);
%     for i=1:length(k_holiday)
%         K3(k_holiday(i))=1;
%     end
% K3=ones(24,1)*K3;K3=K3(:);K3=K3';
% 
% K=[K1;K2;K3];
% K=K(:,(startday-1)*24+1:24*d);

startday=71;
d=104;
data = NTL((startday-1)*96+1:96*d,:)'; %data(1,:)=sum(data,1);
datapv = PV((startday-1)*96+1:96*d,:)';

% data = chickenpox_dataset;data = [data{:}];
 numTimeStepsTrain = floor(0.9*(length(data)));

dataTrain = data(:,1:numTimeStepsTrain+1);
dataTest = data(:,numTimeStepsTrain+1:end);
%k=size(dataTrain);k=min(k);
        mu=zeros(k,1);sig=zeros(k,1);

        for i=1:k
            mu(i) = mean(dataTrain(i,:));
            sig(i) = std(dataTrain(i,:));
        end
        dataTrainStandardized = (dataTrain - mu) ./ sig;
        XTrain = dataTrainStandardized(1:k,1:end-1);%remain1
        YTrain = dataTrainStandardized(1:k,2:end);% emit first one
        numFeatures = k; %
        numResponses = k; %
        numHiddenUnits = 200;

        layers = [ ...
            sequenceInputLayer(numFeatures)
            lstmLayer(numHiddenUnits)
            fullyConnectedLayer(numResponses)
            regressionLayer];
        options = trainingOptions('adam', ...
            'MaxEpochs',250, ...
            'GradientThreshold',1, ...
            'InitialLearnRate',0.005, ...
            'LearnRateSchedule','piecewise', ...
            'LearnRateDropPeriod',125, ...
            'LearnRateDropFactor',0.2, ...
            'Verbose',0, ...
            'Plots','training-progress');
        net = trainNetwork(XTrain,YTrain,layers,options);
%        load lstmnet0213
        dataTestStandardized = (dataTest(1:k,:) - mu(1:k,:)) ./ sig(1:k,:);
        XTest = dataTestStandardized(:,1:end-1);
        net = predictAndUpdateState(net,XTrain);%update last point
        [net,YPred] = predictAndUpdateState(net,YTrain(:,end));

        numTimeStepsTest = numel(XTest(1,:));
        for i = 2:numTimeStepsTest
            [net,YPred(:,i)] = predictAndUpdateState(net,YPred(:,i-1),'ExecutionEnvironment','cpu');
        end
        YPred = sig(1:k,:).*YPred + mu(1:k,:);
        YTest = dataTest(:,2:end);
        rmse1 = sqrt(mean((YPred(1,:)-YTest(1,:)).^2));
        nrmse1=rmse1/(max(YTest(1,:))-min(YTest(1,:)));        
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
% 
        %resetnet
        net = resetState(net);
        net = predictAndUpdateState(net,XTrain);
        YPred = [];
        numTimeStepsTest = numel(XTest(1,:));
        XTest(1,:)=XTest(1,:).*(0.8+0.2*rand(numTimeStepsTest,1))';
        for i = 1:numTimeStepsTest
            [net,YPred(:,i)] = predictAndUpdateState(net,XTest(:,i),'ExecutionEnvironment','cpu');
        end
        YPred = sig(1:k,:).*YPred + mu(1:k,:);
        rmse2 = sqrt(mean((YPred(1,:)-YTest(1,:)).^2));
        nrmse2=rmse2/(max(YTest(1,:))-min(YTest(1,:)));
        Error_ind(k,:)=[rmse2 nrmse2];
%%        
        
        subplot(2,1,2)
        plot(YTest(1,:))
        hold on
        plot(YPred(1,:),'.-')
        hold off
        legend(["Observed" "Predicted"])
        ylabel("Cases")
        title("Forecast with Updates")
        clearvars -except Error_ind
end
%end

