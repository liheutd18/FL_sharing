% LSTM_CLEAN
clear
[ntl l pv] = Pecan();
aggntl = sum(ntl,2);

[col row_idx] = size (ntl);

load('Results0606.mat')   
%data = chickenpox_dataset;
%data = [data{:}];
% figure
% plot(data)
% xlabel("Month")
% ylabel("Cases")
% title("Monthy Cases of Chickenpox")

data = ntl(1:312,:)';

%%%% Data correlation analysis
%Alldata = [data; ntl(1:312,:)'];
YPred_all = zeros(30,71);

for j = 1:10
numTimeStepsTrain = floor(0.77*numel(data(j,:)));

dataTrain = data(j,1:numTimeStepsTrain+1);
dataTest = data(j,numTimeStepsTrain+1:end);

mu = mean(dataTrain);
sig = std(dataTrain);

dataTrainStandardized = (dataTrain - mu) / sig;

XTrain = dataTrainStandardized(1:end-1);
YTrain = dataTrainStandardized(2:end);

% numFeatures = 1;
% numResponses = 1;
% numHiddenUnits = 100;
% 
% layers = [ ...
%     sequenceInputLayer(numFeatures)
%     lstmLayer(numHiddenUnits)
%     fullyConnectedLayer(numResponses)
%     regressionLayer];
% 
% options = trainingOptions('adam', ...
%     'MaxEpochs', 150, ...
%     'GradientThreshold',1, ...
%     'InitialLearnRate',0.005, ...
%     'LearnRateSchedule','piecewise', ...
%     'LearnRateDropPeriod',125, ...
%     'LearnRateDropFactor',0.2, ...
%     'Verbose',0);%,...    'Plots','training-progress');
% 
% net = trainNetwork(XTrain,YTrain,layers,options);

net = networks{j};
net1 = net;


%% Update weights OR NOT

% netFL = weightsupdatewobias(networks,numHiddenUnits);
% net = netFL;

%%%% Forecast Future Time Steps

dataTestStandardized = (dataTest - mu) / sig;
XTest = dataTestStandardized(1:end-1);

%%%% Select XTrain: limited access of historical
%XTrain = XTrain (1,217:240);

%%%% Hidden state changes
net = resetState(net);
%net = predictAndUpdateState(net,XTrain);
[net,YPred] = predictAndUpdateState(net,YTrain(end));

%%%% Check update
%net.Layers(2, 1).InputWeights - net1.Layers(2, 1).InputWeights  

numTimeStepsTest = numel(XTest);
for i = 2:numTimeStepsTest
    [net,YPred(:,i)] = predictAndUpdateState(net,YPred(:,i-1),'ExecutionEnvironment','cpu');
end

YPred = sig*YPred + mu;

YTest = dataTest(2:end);
rmse1 = sqrt(mean((YPred-YTest).^2));

YPred_all((j-1)*3+1,:) = YPred;
YPred_all((j-1)*3+1+2,:) = YTest;

% figure
% plot(dataTrain(1:end-1))
% hold on
% idx = numTimeStepsTrain:(numTimeStepsTrain+numTimeStepsTest);
% plot(idx,[data(numTimeStepsTrain) YPred],'.-')
% hold off
% xlabel("Time")
% ylabel("Power")
% title("Forecast")
% legend(["Observed" "Forecast"])

%%%% figure
% subplot(2,1,1)
% plot(YTest)
% xlim([0 72])
% hold on
% plot(YPred,'.-')
% hold off
% legend(["Observed" "Forecast"])
% ylabel("Power (kW)")
% title("Forecast")



%%%% Update Network State with Observed Values
% 
net = resetState(net);
net = predictAndUpdateState(net,XTrain);

YPred = [];
numTimeStepsTest = numel(XTest);
for i = 1:numTimeStepsTest
    [net,YPred(:,i)] = predictAndUpdateState(net,XTest(:,i),'ExecutionEnvironment','cpu');
end

YPred = sig*YPred + mu;

rmse2 = sqrt(mean((YPred-YTest).^2));
YPred_all((j-1)*3+1+1,:) = YPred;

%%%% figure
% subplot(2,1,2)
% plot(YTest)
% xlim([0 72])
% hold on
% plot(YPred,'.-')
% hold off
% legend(["Observed" "Predicted"])
% ylabel("Power (kW)")
% title("Forecast with Updates")
RMSE_all(j,:) = [rmse1 rmse2];
max_min = max(data')-min(data');

end

% subplot(2,1,2)
% stem(YPred - YTest)
% xlabel("Time (hour)")
% ylabel("Error")
% title("RMSE = " + rmse)