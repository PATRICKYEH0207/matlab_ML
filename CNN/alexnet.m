%alenxet
path=data_path;
figuresdir_AUC=AUC_path;
figuresdir_ConfusionMatrix=onfusionMatrix_path;
figuresdir_erro=erro_path;
%set option
MiniBatchSize=Num(1,1);Epochs=Num(1,2);
InitialLearnRate=Num(1,3);ValidationData=Num(1,4);
ValidationFrequency=Num(1,5);transfer=Num(1,6);
Augmenter=Num(1,7);times=Num(1,8);Verbose=Num(1,9);
Solver=string(txt(1,1));Plots=string(txt(1,2));
environment=string(txt(1,3));Shuffle=string(txt(1,4));
finalcsv={'User Name',path;...
        'Solver',Solver;...
        'environment',environment;...
        'Shuffle',Shuffle;...
        'Verbose',Verbose;...
        'transfer',transfer;...
        'Augmenter',Augmenter;...
        'Epochs',Epochs;'times',times;...
        'MiniBatchSize',MiniBatchSize;...
        'ValidationFrequency',ValidationFrequency;...
        'InitialLearnRate',InitialLearnRate;...
        'ValidationData',ValidationData;...
        };
File = cellstr(path);
imds = imageDatastore(File, ...
    'IncludeSubfolders',true, ...
    'LabelSource','foldernames');
for time=1:times
        [imdsTrain,imdsValidation] = splitEachLabel(imds,ValidationData,'randomized');
        if transfer==1
            net = alexnet;
            inputSize = net.Layers(1).InputSize;
            layers = net.Layers;
            layers(23) = fullyConnectedLayer(2);
            layers(25) = classificationLayer;
        else
               layers = [
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
        end
        pixelRange = [-5 5];
        if Augmenter==1
            imageAugmenter = imageDataAugmenter( ...
                'RandXTranslation',pixelRange, ...
                'RandXReflection',true, ...
                'RandRotation',[0 360], ...
                'RandScale',[0.5 1],...
                'RandYTranslation',pixelRange);            
            augimdsTrain = augmentedImageDatastore(inputSize(1:2),imdsTrain, ...
                        'DataAugmentation',imageAugmenter);
            augimdsValidation = augmentedImageDatastore(inputSize(1:2),imdsValidation);
        else
            augimdsTrain = augmentedImageDatastore(inputSize(1:2),imdsTrain);
            augimdsValidation = augmentedImageDatastore(inputSize(1:2),imdsValidation);
        end
        %-----traning-----%
        Verbose=logical(mod(Verbose,0));
        options = trainingOptions(Solver, ...
                'MiniBatchSize',MiniBatchSize, ...
                'MaxEpochs',Epochs, ...
                'InitialLearnRate',InitialLearnRate, ...
                'Shuffle',Shuffle, ...
                'ValidationData',augimdsValidation, ...
                'ValidationFrequency',ValidationFrequency, ...
                'Verbose',Verbose, ...
                'ExecutionEnvironment',environment,...
                'OutputFcn',@(info)savetrainingplot(info),...
                'Plots',Plots);%
            diary DiaryFile.txt
            netTransfer = trainNetwork(augimdsTrain,layers,options);
            %savefig(netTransfer,'training-progress.png')
            diary off
            %type myDiaryFile.xlsx
            [YPred,scores] = classify(netTransfer,augimdsValidation);
            YValidation = imdsValidation.Labels;
            finalcsv(16,time+1) = num2cell(mean(YPred == YValidation));
    %erro image output
            for i=1:augimdsValidation.NumObservations
                file_erro = char(imdsValidation.Files(i));
                [filepath_erro,name,ext] = fileparts(file_erro);
                test_Name=cellstr(YPred(i));Validation_Name=cellstr(YValidation(i));
                Compare=strcmp(test_Name,Validation_Name);
                if Compare==0
                    image=imread(file_erro);
                    if transfer==1
                        imwrite(image,['C:\xampp\htdocs\transfer\erro\Alex_',num2str(time),'_',char(Validation_Name),'_',name,ext]);
                    else
                        imwrite(image,['C:\xampp\htdocs\transfer\erro\nonAlex_',num2str(time),'_',char(Validation_Name),'_',name,ext]);
                    end
                end
            end
            %-----confusion matrix-----%
            C = confusionmat(YValidation,YPred);
            plotconfusion(YValidation,YPred)
            if transfer==1
                filename=['ConfusionMatrix_Alex_',num2str(time),'.png'];
                saveas(gcf,strcat(figuresdir_ConfusionMatrix,filename))
            else
                filename=['ConfusionMatrix_nonAlex_',num2str(time),'.png'];
                saveas(gcf,strcat(figuresdir_ConfusionMatrix,filename))
            end
            TP=C(1,1);FN=C(2,1);FP=C(1,2);TN=C(2,2);
            TPR=TP/(TP+FN);
            FPR=FP/(FP+TN);
            finalcsv{14,1}='sensitivity';
            finalcsv{15,1}='1-specificity';
            finalcsv(14,time+1)=num2cell(TPR);
            finalcsv(15,time+1)=num2cell(TPR);
%-----roc cruve-----% 
            YPred = predict(netTransfer, augimdsValidation);
            YTest = zeros(size(YPred));
            scores = zeros(size(imdsValidation.Labels));
            labels = zeros(size(imdsValidation.Labels));
            a = dir(char(File));
            for n=1:length(imdsValidation.Labels)
                scores(n) = YPred(n,1);
                    if imdsValidation.Labels(n) == a(length(a)-1).name
                        YTest(n,1) = 1;
                        labels(n) = true;
                    elseif imdsValidation.Labels(n) == a(length(a)).name
                        YTest(n,2) = 1;
                        labels(n) = false;
                    end
            end
            [X,Y,T,AUC] = perfcurve(labels==1,scores,'true');
            figure
            plot(X,Y)
            xlabel('False positive rate') 
            ylabel('True positive rate')
            if transfer==1
                filename=['AUC_Alex_',num2str(time),'.png'];
                saveas(gcf,strcat(figuresdir_AUC,filename))
            end
            if transfer==0
                filename=['AUC_nonAlxe_',num2str(time),'.png'];
                saveas(gcf,strcat(figuresdir_AUC,filename))
            end
end
finalcsv{16,1}='ACC';
%creat csv
if transfer==1
    fid = fopen('C:\xampp\htdocs\transfer\ACC\alexnet.csv','w');
    for i=1:4
        fprintf(fid,'%s,%s\n',finalcsv{i,1},finalcsv{i,2});
    end
    for i=5:11
        fprintf(fid,'%s,%i\n',finalcsv{i,1},finalcsv{i,2});
    end
	for i=12:13
        fprintf(fid,'%s,%f\n',finalcsv{i,1},finalcsv{i,2});
    end
	for i=14:16
        fprintf(fid,'%s,',finalcsv{i,1});
        for j=2:time+1
            fprintf(fid,'%f,',finalcsv{i,j});
        end
        fprintf(fid,'\n');
	end
    fclose(fid);
    close all;
else 
    fid = fopen('nonalexnet.csv','w');
     for i=1:5
        fprintf(fid,'%s,%s\n',finalcsv{i,1},finalcsv{i,2});
    end
    for i=6:11
        fprintf(fid,'%s,%i\n',finalcsv{i,1},finalcsv{i,2});
    end
	for i=12:13
        fprintf(fid,'%s,%f\n',finalcsv{i,1},finalcsv{i,2});
    end
	for i=14:16
        fprintf(fid,'%s,',finalcsv{i,1});
        for j=2:time+1
            fprintf(fid,'%f,',finalcsv{i,j});
        end
        fprintf(fid,'\n');
    end
    fclose(fid); 
end
