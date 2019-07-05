function [inputSize,lgraph]=net(net_name,imdsTrain)
    if net_name==1
        net = alexnet;
        inputSize = net.Layers(1).InputSize;
        lgraph = net.Layers;
        lgraph(23) = fullyConnectedLayer(2);
        lgraph(25) = classificationLayer;
    elseif net_name==2
        lgraph = [
                imageInputLayer([227 227 3])%1

                convolution2dLayer(11,96,'Stride',4)%2
                reluLayer%3
                crossChannelNormalizationLayer(5)%4
                maxPooling2dLayer(3,'Stride',2)%5

                convolution2dLayer(5,256,'Stride',1,'Padding',2)%6
                reluLayer%7
                crossChannelNormalizationLayer(5)%8
                maxPooling2dLayer(3,'Stride',2)%9

                convolution2dLayer(3,384,'Stride',1,'Padding',1)%10
                reluLayer%11
                convolution2dLayer(3,384,'Stride',1,'Padding',2)%12
                reluLayer%13
                convolution2dLayer(3,256,'Stride',1,'Padding',2)%14
                reluLayer%15
                maxPooling2dLayer(3,'Stride',2)%16
                fullyConnectedLayer(4096)%17
                reluLayer%18
                dropoutLayer%19
                fullyConnectedLayer(4096)%20
                reluLayer%21
                dropoutLayer%22
                fullyConnectedLayer(2)%23
                softmaxLayer%24
                classificationLayer%25
        ];
        inputSize = layers(1,1).InputSize;
    elseif net_name==3
        net = inceptionv3;
        inputSize = net.Layers(1).InputSize;
        % Extract the layer graph from the trained network and plot the layer graph.
        lgraph = layerGraph(net);
        % Replace Final Layers
        lgraph = removeLayers(lgraph, {'predictions','predictions_softmax','ClassificationLayer_predictions'});
        %修改最後fc,sofmax,Classification名稱
        numClasses = numel(categories(imdsTrain.Labels));
        newLayers = [
        fullyConnectedLayer(numClasses,'Name','fc','WeightLearnRateFactor',10,'BiasLearnRateFactor',10)
        softmaxLayer('Name','softmax')
        classificationLayer('Name','classoutput')];
        lgraph = addLayers(lgraph,newLayers);
        lgraph = connectLayers(lgraph,'avg_pool','fc');
    elseif net_name==4
        net = resnet101;
        lgraph = layerGraph(net);
        lgraph = removeLayers(lgraph, {'fc1000','prob','ClassificationLayer_predictions'});
        numClasses = numel(categories(imdsTrain.Labels));
        newLayers = [
            fullyConnectedLayer(numClasses,'Name','fc','WeightLearnRateFactor',10,'BiasLearnRateFactor',10)
            softmaxLayer('Name','softmax')
            classificationLayer('Name','classoutput')];
        lgraph = addLayers(lgraph,newLayers);
        lgraph = connectLayers(lgraph,'pool5','fc');
    end
end