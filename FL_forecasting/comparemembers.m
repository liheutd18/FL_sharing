clear

%%

load('Results0606.mat')

for combination = 1:5

for members = 1:10
netFL = weightsupdatewobiasmembers(networks,numHiddenUnits,members,combination);
net = netFL;

%%%% Forecast Future Time Steps

dataTestStandardized = (dataTest - mu) / sig;
XTest = dataTestStandardized(1:end-1);

%%%% Select XTrain: limited access of historical
%XTrain = XTrain (1,217:end);

%%%% Hidden state changes
net = resetState(net);
net = predictAndUpdateState(net,XTrain);
[net,YPred] = predictAndUpdateState(net,YTrain(end));

numTimeStepsTest = numel(XTest);
for i = 2:numTimeStepsTest
    [net,YPred(:,i)] = predictAndUpdateState(net,YPred(:,i-1),'ExecutionEnvironment','cpu');
end

YPred = sig*YPred + mu;

YTest = dataTest(2:end);
rmse1(combination,members) = sqrt(mean((YPred-YTest).^2));

net = resetState(net);
net = predictAndUpdateState(net,XTrain);

YPred = [];
numTimeStepsTest = numel(XTest);
for i = 1:numTimeStepsTest
    [net,YPred(:,i)] = predictAndUpdateState(net,XTest(:,i),'ExecutionEnvironment','cpu');
end

YPred = sig*YPred + mu;

rmse2(combination,members) = sqrt(mean((YPred-YTest).^2));


end
end


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


% plot(YTest)
% hold on
% plot(YPred,'.-')
% hold off
% legend(["Observed" "Forecast"])
% ylabel("Power (kW)")
% title("Forecast")

% subplot(2,1,2)
% stem(YPred - YTest)
% xlabel("Time (hour)")
% ylabel("Error")
% title("RMSE = " + rmse)

%%%% Update Network State with Observed Values
  
% net = resetState(net);
% net = predictAndUpdateState(net,XTrain);
% 
% YPred = [];
% numTimeStepsTest = numel(XTest);
% for i = 1:numTimeStepsTest
%     [net,YPred(:,i)] = predictAndUpdateState(net,XTest(:,i),'ExecutionEnvironment','cpu');
% end
% 
% YPred = sig*YPred + mu;
% 
% rmse2 = sqrt(mean((YPred-YTest).^2));
% 
% %figure
% subplot(4,1,2)
% plot(YTest)
% hold on
% plot(YPred,'.-')
% hold off
% legend(["Observed" "Predicted"])
% ylabel("Power (kW)")
% title("Forecast with Updates")
% 
% load('Results0606onenet.mat')
% 
% dataTestStandardized = (dataTest - mu) / sig;
% XTest = dataTestStandardized(1:end-1);
% 
% %%%% Select XTrain: limited access of historical
% %XTrain = XTrain (1,217:end);
% 
% %%%% Hidden state changes
% net = resetState(net);
% net = predictAndUpdateState(net,XTrain);
% [net,YPred] = predictAndUpdateState(net,YTrain(end));
% 
% numTimeStepsTest = numel(XTest);
% for i = 2:numTimeStepsTest
%     [net,YPred(:,i)] = predictAndUpdateState(net,YPred(:,i-1),'ExecutionEnvironment','cpu');
% end
% 
% YPred = sig*YPred + mu;
% 
% YTest = dataTest(2:end);
% rmse3 = sqrt(mean((YPred-YTest).^2));
% 
% % figure
% % plot(dataTrain(1:end-1))
% % hold on
% % idx = numTimeStepsTrain:(numTimeStepsTrain+numTimeStepsTest);
% % plot(idx,[data(numTimeStepsTrain) YPred],'.-')
% % hold off
% % xlabel("Time")
% % ylabel("Power")
% % title("Forecast")
% % legend(["Observed" "Forecast"])
% 
% %figure
% subplot(4,1,3)
% plot(YTest)
% hold on
% plot(YPred,'.-')
% hold off
% legend(["Observed" "Forecast"])
% ylabel("Power (kW)")
% title("Forecast")
% 
% % subplot(2,1,2)
% % stem(YPred - YTest)
% % xlabel("Time (hour)")
% % ylabel("Error")
% % title("RMSE = " + rmse)
% 
% %%%% Update Network State with Observed Values
% % 
% net = resetState(net);
% net = predictAndUpdateState(net,XTrain);
% 
% YPred = [];
% numTimeStepsTest = numel(XTest);
% for i = 1:numTimeStepsTest
%     [net,YPred(:,i)] = predictAndUpdateState(net,XTest(:,i),'ExecutionEnvironment','cpu');
% end
% 
% YPred = sig*YPred + mu;
% 
% rmse4 = sqrt(mean((YPred-YTest).^2));
% 
% %figure
% subplot(4,1,4)
% plot(YTest)
% hold on
% plot(YPred,'.-')
% hold off
% legend(["Observed" "Predicted"])
% ylabel("Power (kW)")
% title("Forecast with Updates")
% 
% RMSE_ALL = [rmse1 rmse2 rmse3 rmse4];