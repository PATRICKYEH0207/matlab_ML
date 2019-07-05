clear all;close all;clc;
[FileName,PathName,FilterIndex] = uigetfile({'*.jpg'},'Select file(s)','MultiSelect','on');
Size = size(FileName);

for x = 1:Size(2)
    File = strcat(PathName, FileName(x));
    %DoShapeTexture(File, FileName(x),2);
    %File=string(File);
    gaborArray = gaborFilterBank(5,8,39,39);  
    %featureVector = gaborFeatures(File, FileName(x),gaborArray,4,4,2);
end
