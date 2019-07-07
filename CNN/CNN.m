clc;clear;
path='';
figuresdir_AUC='';
figuresdir_ConfusionMatrix='';
filename = 'option_num.csv';
Num = readmatrix(filename);
MiniBatchSize=Num(1,1);Epochs=Num(1,2);
InitialLearnRate=Num(1,3);ValidationData=Num(1,4);
ValidationFrequency=Num(1,5);net_name=Num(1,6);
Augmenter=Num(1,7);times=Num(1,8);Verbose=Num(1,9);
filename = 'option_str.csv';
[num,txt,raw] =xlsread(filename);
Solver=string(txt(1,1));Plots=string(txt(1,2));
environment=string(txt(1,3));Shuffle=string(txt(1,4));
finalcsv={'User Name',path;...
        'Solver',Solver;...
        'environment',environment;...
        'Shuffle',Shuffle;...
        'Verbose',Verbose;...
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
imds.ReadFcn = @(imgname) imresize(imread(imgname), [299,299]);
for time=1:times
        [imdsTrain,imdsValidation] = splitEachLabel(imds,ValidationData,'randomized');
        [inputSize,lgraph]=net(net_name,imdsTrain);
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
            netTransfer = trainNetwork(augimdsTrain,lgraph,options);
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
                    imwrite(image,['C:\...',num2str(time),'_',char(Validation_Name),'_',name,ext]);
                end
            end
            %-----confusion matrix-----%
            C = confusionmat(YValidation,YPred);
            plotconfusion(YValidation,YPred)
            filename=['ConfusionMatrix_',num2str(time),'.png'];
            saveas(gcf,strcat(figuresdir_ConfusionMatrix,filename))
            TP=C(1,1);FN=C(2,1);FP=C(1,2);TN=C(2,2);
            TPR=TP/(TP+FN);%sensitivity±Ó·P«×
            FPR=FP/(FP+TN);%1-specificity¯S²§«×
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
            filename=['AUC_',num2str(time),'.png'];
            saveas(gcf,strcat(figuresdir_AUC,filename))
end
finalcsv{16,1}='ACC';
%creat csv
fid = fopen('C:\...\net.csv','w');
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
