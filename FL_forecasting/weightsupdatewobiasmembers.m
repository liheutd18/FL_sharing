function net_new = weightsupdatewobiasmembers(networks,numHiddenUnits,members,combination)

%%%% Learning rate eta
eta_FL = 0.001;
%%%% Choose target networks
if combination == 1 
networks = { networks{members}};
elseif combination == 2 
networks = { networks{6} networks{members}  };
elseif combination == 3 
networks = { networks{6} networks{2} networks{members} };
elseif combination == 4 
networks = { networks{6} networks{2} networks{1} networks{members}};
elseif combination == 5 
networks = { networks{6} networks{2} networks{1} networks{7} networks{members}};
end
%networks = {networks{1} networks{2} networks{3} networks{4} networks{6} networks{7} networks{9}};

for i = 1 : length(networks)
    
%%%% networks and net are the local net
    net = networks{i};
    
%%%% Local update weights
    layers_new = net.Layers;
    
    layers_new(2,1).InputWeights = layers_new(2,1).InputWeights - eta_FL.*gradient(net.Layers(2, 1).InputWeights);
    layers_new(2,1).RecurrentWeights = layers_new(2,1).RecurrentWeights - eta_FL.*gradient(net.Layers(2, 1).RecurrentWeights);
    layers_new(2,1).Bias = layers_new(2,1).Bias - eta_FL.*gradient(net.Layers(2, 1).Bias);

    weights_agg = [layers_new(2,1).InputWeights, layers_new(2,1).RecurrentWeights, layers_new(2,1).Bias]; % 400*1 + 400*100 + 400*1
    Weights_table{i} = weights_agg;

end

%%
%%%% Global update weights
%%%% layers_new_all is the global net
    % just copy a net format
    layers_new_all = net.Layers;
    
    % initilize the net values to 0 or random?     
    layers_new_all(2,1).InputWeights = single(zeros(numHiddenUnits*4,1));
    layers_new_all(2,1).RecurrentWeights = single(zeros(numHiddenUnits*4,numHiddenUnits));
    weights_agg_new_all{1} = [layers_new_all(2,1).InputWeights, layers_new_all(2,1).RecurrentWeights];
    
   
%%%% interation of j_max rounds
    iter_max = 50;
    error = zeros(1,iter_max);
    
for iter = 2:iter_max    
    
    for i = 1:length(networks)
        layers_new_all(2,1).InputWeights = layers_new_all(2,1).InputWeights + Weights_table{i}(:,1);
        layers_new_all(2,1).RecurrentWeights = layers_new_all(2,1).RecurrentWeights + Weights_table{i}(:,2:numHiddenUnits+1);
    end
    
    layers_new_all(2,1).InputWeights = layers_new_all(2,1).InputWeights/length(networks);
    layers_new_all(2,1).RecurrentWeights = layers_new_all(2,1).RecurrentWeights/length(networks);
    
    
    weights_agg_new_all{iter} = [layers_new_all(2,1).InputWeights, layers_new_all(2,1).RecurrentWeights];
    error(iter) = sum(sum(abs(weights_agg_new_all{iter} - weights_agg_new_all{iter-1})));
    
    if error(iter) < 1E-3
       break 
    end
end  

%%%% Modification made in net
    modify_able_net = net.saveobj;
    modify_able_net.Layers = layers_new_all; %layers_new_all
    Modified_net = net.loadobj(modify_able_net);
    net_new = Modified_net; % net_new.Layers(2, 1).InputWeights  - net.Layers(2, 1).InputWeights


end

%%%% retrain the network
%net_new = trainNetwork(XTrain,YTrain,layers_new,options);
