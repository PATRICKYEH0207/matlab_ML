[FileName,PathName,FilterIndex] = uigetfile({'*.jpg'},'Select file(s)','MultiSelect','on');
Size = size(FileName);

for x = 1:Size(2)
    File = strcat(PathName, FileName(x));
    %GLCM feature
    DoShapeTexture(File, FileName(x),dimension);
    File=string(File);
    gaborArray = gaborFilterBank(u,v,m,n);  
    featureVector = gaborFeatures(File, FileName(x),gaborArray,d1,d2,dimension);
end
