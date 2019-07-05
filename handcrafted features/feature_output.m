clc; clear; close;
imds = imageDatastore(Filepath, ...%choose data
    'IncludeSubfolders',true, ...
    'LabelSource','foldernames');
net = netname;%choose net
analyzeNetwork(net)
inputSize = net.Layers(1).InputSize;
augimds = augmentedImageDatastore(inputSize(1:2),imds);
layer = 'predictions';
features = activations(net,augimds,layer,'OutputAs','rows');
csvwrite(feature_name,features)%Enter the file name