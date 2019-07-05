%load ionosphere
X=csvread(feature_file,1);
[num, Y] = xlsread(class_file);
ionosphere = array2table(X);
ionosphere.Group = Y;