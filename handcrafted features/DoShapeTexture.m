function [Result] = DoShapeTexture(File, FileName,dimension)
fprintf('%s\n', FileName{1,1});
ImgFile = imread(File{1,1});
ImgFile = ImgFile(:,:,dimension);
%ImgFile = imread('E:\TMU-2Research\Brain\TumorMRI\3-Delineation\GBM\ROI\190955-5-95clear.bmp');
Size = size(ImgFile);
BW = ImgFile;
Area = 0;
for x = 1:Size(1)
    for y = 1:Size(2)
        if (ImgFile(x,y) == 0)
            BW(x,y) = 0;
        else
            BW(x,y) = 255;
            Area = Area+1;
        end
    end
end
%SSum=0;
%Vsum=0;
%S=test(:, :, 2);
%V=test(:, :, 3);
%for x = 1:Size(1)
%    for y = 1:Size(2)
%        SSum=SSum+S;
%       Vsum=Vsum+V;
%    end
%end

%fprintf('%0.2f %0.2f %s\n',SSum, Vsum ,FileName)

%test= imread(File{1,1});
%subplot(3,3,7);
%imshow(test);
%title('original');

%test= double(test)/255;
%r=test(:,:,1);
%g=test(:,:,2);
%b=test(:,:,3);
%k=(r-1).^2+(g-1).^2+(b-1).^2;
%k=exp(-k);
%subplot(3,1,2)
%imshow(k,[]);
%title('test');

%test=rgb2hsv(test);
%H=test(:, :, 1);
%S=test(:, :, 2);
%V=test(:, :, 3);
%subplot(3,3,1);
%imshow(H);
%title('H');
%subplot(3,3,2);
%imshow(S);
%title('S');
%subplot(3,3,3);
%imshow(V);
%title('V');
%title('hsv');
%i=0;
%for x = 1:Size(1)
%    for y = 1:Size(2)
%        if ( V(x,y)>0.95 && S(x,y)<0.45 )
%            S(x,y)=0;
%            V(x,y)=0;
%        end
%    end
%end
%subplot(3,3,4);
%imshow(H);
%title('H');
%subplot(3,3,5);
%imshow(S);
%title('S');
%subplot(3,3,6);
%imshow(V);
%title('V');
%for x = 1:Size(1)
%    for y = 1:Size(2)
%        if (ImgFile(x,y) == 255)
%            V(x,y)=0;
        %test(x,y) =257 - ImgFile(x,y);
%        end
%    end
%end
%subplot(3,1,3);
%imshow(V);
%title('V=255');

BW_filled = imfill(BW,'holes');%對圖像進行填充
Boundary = bwboundaries(BW_filled);%我們將影像轉為黑白影像，以利接下來用bwboundaries來判定邊界。
BoundarySize = size(Boundary{1});
Perimeter = BoundarySize(1);
Compactness = 1-4*3.14159*Area/(Perimeter*Perimeter);%緊密度(compactness)：是屬於形狀上的特徵值(shape features)，藉由影像物體的周長與面積的關係，來描其物體的圓形曲率與緊密度的外形特徵。其定義為周長平方與面積之比率，當物體愈圓，則比率接近1。 
%定義表示法:：C = (周長平方)/4π(面積),C表示緊密度，當物體為圓形時 C 相當於 (2πr)^2/4π(πr^2) = 1。
%NRL mean, NRL SD
Coor = Boundary{1};
MeanX = mean(Coor(:,1));
MeanY = mean(Coor(:,2));
MaxDis = 0;
for x = 1:Perimeter
    Dis = power(power((Coor(x,1)-MeanX),2)+power((Coor(x,2)-MeanY),2),0.5);
    if (Dis > MaxDis)
        MaxDis = Dis;
    end
end
NRL_Avg = 0;
for x = 1:Perimeter
    Dis = power(power((Coor(x,1)-MeanX),2)+power((Coor(x,2)-MeanY),2),0.5);
    NRL_Avg = NRL_Avg + double(Dis/MaxDis);
end    
NRL_Avg = NRL_Avg/Perimeter;
NRL_SD = 0;
for x = 1:Perimeter
    Dis = power(power((Coor(x,1)-MeanX),2)+power((Coor(x,2)-MeanY),2),0.5);
    NRL_SD = NRL_SD + power(Dis/MaxDis-NRL_Avg,2);
end    
NRL_SD = power(NRL_SD/Perimeter,0.5);
% for k=1
%    b = boundaries{k};
%    plot(b(:,2),b(:,1),'g','LineWidth',3);
% end
%Name = strcat(File{1,1}, '_Area.bmp');
%imwrite(uint8(BW_filled), Name);
%imshow(BW)
%Record tumor area
Mask = ImgFile;%whatever value works.
%mask = ImgFile >300;
%fixedImage = regionfill(ImgFile, mask);
%subplot(2,1,1)
%imshow(fixedImage);
%subplot(2,1,2)
Max = 0; Min = 255;
MaskNum = 0;
for x = 1:Size(1)
    for y = 1:Size(2)
        if (ImgFile(x,y) == 0)
            Mask(x,y) =  NaN(1);
        else
            Mask(x,y) = ImgFile(x,y);
            MaskNum = MaskNum+1;
            if (Mask(x,y) > Max)
                Max = Mask(x,y);
            end
            if (Mask(x,y) < Min)
                Min = Mask(x,y);
            end
        end
    end
end

             
% Stretch
% for x = 1:Size(1)
%     for y = 1:Size(2)
%         if (ImgFile(x,y) ~= 0)
%             Mask(x,y) = 255*(double(Mask(x,y)-Min)/double(Max-Min));
%         end
%     end
% end
%imshow(Mask);
%LBP
 nFiltSize=8;
% nFiltRadius=1;
% filtR=generateRadialFilterLBP(nFiltSize, nFiltRadius);
% effLBP= efficientLBP(ImgFile, 'filtR', filtR, 'isRotInv', false, 'isChanWiseRot', false);
% %imshow(effLBP);
% for x = 1:Size(1)
%     for y = 1:Size(2)
%         if (ImgFile(x,y) == 0)
%             Mask(x,y) = NaN(1);
%         else
%             Mask(x,y) = effLBP(x,y);
%         end
%     end
% end
%imwrite(uint8(Mask),'E:\TMU-2Research\_Journal_\105_GBM-3_LBP+GLCM\Figures\Original\LBP.bmp');
%imshow(Mask);
% %Histogram
% MaskHis = zeros(MaskNum,1);
% k = 1;
% for x = 1:Size(1)
%     for y = 1:Size(2)
%         if(Mask(x,y) ~= 0)
%             MaskHis(k,:)=Mask(x,y);
%             k = k+1;
%         end
%     end
% end
binsRange=(1:2^nFiltSize)-1;
MaskHis = hist(single( Mask(:) ), binsRange);
MaskHis(:,1)=0;
MeanGray=mean(MaskHis);
VarGray=var(MaskHis);
SkewnessGray=skewness(MaskHis);
KurtosisGray=kurtosis(MaskHis);
%GLCM
kernel = [0 1; 0 -1; 1 0; -1 0]; %[0 1; 0 -1; 1 0; -1 0]
glcmOri = graycomatrix( Mask, 'NumLevels', 8, 'GrayLimits', [],'offset', kernel );
FeaturesOri = GLCM_Features1(glcmOri,0);
FFile = fopen('.\GLCMB_originalG.csv', 'a');
fprintf(FFile, '%s, ', FileName{1,1});
fprintf(FFile, '%2.5f, ', FeaturesOri.autoc(1));
fprintf(FFile, '%2.5f, ', FeaturesOri.contr(1));
fprintf(FFile, '%2.5f, ', FeaturesOri.corrm(1));
fprintf(FFile, '%2.5f, ', FeaturesOri.corrp(1));
fprintf(FFile, '%2.5f, ', FeaturesOri.cprom(1));
fprintf(FFile, '%2.5f, ', FeaturesOri.cshad(1));
fprintf(FFile, '%2.5f, ', FeaturesOri.dissi(1));
fprintf(FFile, '%2.5f, ', FeaturesOri.energ(1));
fprintf(FFile, '%2.5f, ', FeaturesOri.entro(1));
fprintf(FFile, '%2.5f, ', FeaturesOri.homom(1));
fprintf(FFile, '%2.5f, ', FeaturesOri.homop(1));
fprintf(FFile, '%2.5f, ', FeaturesOri.maxpr(1));
fprintf(FFile, '%2.5f, ', FeaturesOri.sosvh(1));
fprintf(FFile, '%2.5f, ', FeaturesOri.savgh(1));
fprintf(FFile, '%2.5f, ', FeaturesOri.svarh(1));
fprintf(FFile, '%2.5f, ', FeaturesOri.senth(1));
fprintf(FFile, '%2.5f, ', FeaturesOri.dvarh(1));
fprintf(FFile, '%2.5f, ', FeaturesOri.denth(1));
fprintf(FFile, '%2.5f, ', FeaturesOri.inf1h(1));
fprintf(FFile, '%2.5f, ', FeaturesOri.inf2h(1));
fprintf(FFile, '%2.5f, ', FeaturesOri.indnc(1));
fprintf(FFile, '%2.5f, ', FeaturesOri.idmnc(1));
fprintf(FFile, '%2.5f, ', MeanGray);
fprintf(FFile, '%2.5f, ', VarGray);
fprintf(FFile, '%2.5f, ', SkewnessGray);
fprintf(FFile, '%2.5f, ', KurtosisGray);
fprintf(FFile, '%2.5f, ', Area);
fprintf(FFile, '%2.5f, ', Perimeter);
fprintf(FFile, '%2.5f, ', Compactness);
fprintf(FFile, '%2.5f, ', NRL_Avg);
fprintf(FFile, '%2.5f, ', NRL_SD);
fprintf(FFile, '\n');
fclose(FFile);
%Ranklet GLCM
% [vertical horizontal diagonal] = RankletFilter(File{1,1},4,4)
% vertical = abs(vertical)*255;
% horizontal = abs(horizontal)*255;
% diagonal = abs(diagonal)*255;
% %imwrite(uint8(vertical),'E:\研究\Brain\TumorMRI\Ranklet_Matlab\RV.bmp');
% %imwrite(uint8(horizontal),'E:\研究\Brain\TumorMRI\Ranklet_Matlab\RH.bmp');
% %imwrite(uint8(diagonal),'E:\研究\Brain\TumorMRI\Ranklet_Matlab\RD.bmp');
% for x = 1:Size(1)
%     for y = 1:Size(2)
%         if Mask(x,y) == NaN;
%             vertical(x,y) = NaN;
%             horizontal(x,y) = NaN;
%             diagonal(x,y) = NaN;
%         end
%     end
% end
% glcmV = graycomatrix( vertical, 'NumLevels', 8, 'GrayLimits', [],'offset', kernel );
% glcmH = graycomatrix( horizontal, 'NumLevels', 8, 'GrayLimits', [],'offset', kernel );
% glcmD = graycomatrix( diagonal, 'NumLevels', 8, 'GrayLimits', [],'offset', kernel );
% %stats = graycoprops( glcm, 'Contrast Correlation Energy Homogeneity');
% FeaturesV = GLCM_Features1(glcmV,0);
% FeaturesH = GLCM_Features1(glcmH,0);
% FeaturesD = GLCM_Features1(glcmD,0);
% fprintf(FFile, '%2.5f, ', FeaturesV.autoc(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.contr(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.corrm(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.corrp(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.cprom(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.cshad(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.dissi(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.energ(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.entro(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.homom(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.homop(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.maxpr(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.sosvh(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.savgh(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.svarh(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.senth(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.dvarh(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.denth(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.inf1h(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.inf2h(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.indnc(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.idmnc(1));
% 
% fprintf(FFile, '%2.5f, ', FeaturesH.autoc(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.contr(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.corrm(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.corrp(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.cprom(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.cshad(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.dissi(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.energ(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.entro(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.homom(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.homop(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.maxpr(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.sosvh(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.savgh(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.svarh(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.senth(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.dvarh(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.denth(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.inf1h(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.inf2h(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.indnc(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.idmnc(1));
% 
% fprintf(FFile, '%2.5f, ', FeaturesD.autoc(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.contr(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.corrm(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.corrp(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.cprom(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.cshad(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.dissi(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.energ(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.entro(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.homom(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.homop(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.maxpr(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.sosvh(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.savgh(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.svarh(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.senth(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.dvarh(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.denth(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.inf1h(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.inf2h(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.indnc(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.idmnc(1));
% 
% [vertical horizontal diagonal] = RankletFilter(File{1,1},8,8)
% vertical = abs(vertical)*255;
% horizontal = abs(horizontal)*255;
% diagonal = abs(diagonal)*255;
% %imwrite(uint8(vertical),'E:\研究\Brain\TumorMRI\Ranklet_Matlab\RV.bmp');
% %imwrite(uint8(horizontal),'E:\研究\Brain\TumorMRI\Ranklet_Matlab\RH.bmp');
% %imwrite(uint8(diagonal),'E:\研究\Brain\TumorMRI\Ranklet_Matlab\RD.bmp');
% for x = 1:Size(1)
%     for y = 1:Size(2)
%         if Mask(x,y) == NaN;
%             vertical(x,y) = NaN;
%             horizontal(x,y) = NaN;
%             diagonal(x,y) = NaN;
%         end
%     end
% end
% glcmV = graycomatrix( vertical, 'NumLevels', 8, 'GrayLimits', [],'offset', kernel );
% glcmH = graycomatrix( horizontal, 'NumLevels', 8, 'GrayLimits', [],'offset', kernel );
% glcmD = graycomatrix( diagonal, 'NumLevels', 8, 'GrayLimits', [],'offset', kernel );
% %stats = graycoprops( glcm, 'Contrast Correlation Energy Homogeneity');
% FeaturesV = GLCM_Features1(glcmV,0);
% FeaturesH = GLCM_Features1(glcmH,0);
% FeaturesD = GLCM_Features1(glcmD,0);
% fprintf(FFile, '%2.5f, ', FeaturesV.autoc(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.contr(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.corrm(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.corrp(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.cprom(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.cshad(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.dissi(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.energ(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.entro(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.homom(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.homop(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.maxpr(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.sosvh(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.savgh(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.svarh(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.senth(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.dvarh(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.denth(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.inf1h(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.inf2h(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.indnc(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.idmnc(1));
% 
% fprintf(FFile, '%2.5f, ', FeaturesH.autoc(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.contr(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.corrm(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.corrp(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.cprom(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.cshad(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.dissi(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.energ(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.entro(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.homom(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.homop(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.maxpr(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.sosvh(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.savgh(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.svarh(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.senth(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.dvarh(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.denth(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.inf1h(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.inf2h(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.indnc(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.idmnc(1));
% 
% fprintf(FFile, '%2.5f, ', FeaturesD.autoc(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.contr(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.corrm(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.corrp(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.cprom(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.cshad(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.dissi(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.energ(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.entro(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.homom(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.homop(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.maxpr(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.sosvh(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.savgh(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.svarh(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.senth(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.dvarh(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.denth(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.inf1h(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.inf2h(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.indnc(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.idmnc(1));
% 
% [vertical horizontal diagonal] = RankletFilter(File{1,1},16,16)
% vertical = abs(vertical)*255;
% horizontal = abs(horizontal)*255;
% diagonal = abs(diagonal)*255;
% %imwrite(uint8(vertical),'E:\研究\Brain\TumorMRI\Ranklet_Matlab\RV.bmp');
% %imwrite(uint8(horizontal),'E:\研究\Brain\TumorMRI\Ranklet_Matlab\RH.bmp');
% %imwrite(uint8(diagonal),'E:\研究\Brain\TumorMRI\Ranklet_Matlab\RD.bmp');
% for x = 1:Size(1)
%     for y = 1:Size(2)
%         if Mask(x,y) == NaN;
%             vertical(x,y) = NaN;
%             horizontal(x,y) = NaN;
%             diagonal(x,y) = NaN;
%         end
%     end
% end
% glcmV = graycomatrix( vertical, 'NumLevels', 8, 'GrayLimits', [],'offset', kernel );
% glcmH = graycomatrix( horizontal, 'NumLevels', 8, 'GrayLimits', [],'offset', kernel );
% glcmD = graycomatrix( diagonal, 'NumLevels', 8, 'GrayLimits', [],'offset', kernel );
% %stats = graycoprops( glcm, 'Contrast Correlation Energy Homogeneity');
% FeaturesV = GLCM_Features1(glcmV,0);
% FeaturesH = GLCM_Features1(glcmH,0);
% FeaturesD = GLCM_Features1(glcmD,0);
% fprintf(FFile, '%2.5f, ', FeaturesV.autoc(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.contr(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.corrm(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.corrp(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.cprom(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.cshad(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.dissi(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.energ(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.entro(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.homom(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.homop(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.maxpr(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.sosvh(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.savgh(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.svarh(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.senth(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.dvarh(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.denth(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.inf1h(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.inf2h(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.indnc(1));
% fprintf(FFile, '%2.5f, ', FeaturesV.idmnc(1));
% 
% fprintf(FFile, '%2.5f, ', FeaturesH.autoc(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.contr(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.corrm(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.corrp(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.cprom(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.cshad(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.dissi(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.energ(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.entro(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.homom(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.homop(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.maxpr(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.sosvh(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.savgh(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.svarh(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.senth(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.dvarh(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.denth(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.inf1h(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.inf2h(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.indnc(1));
% fprintf(FFile, '%2.5f, ', FeaturesH.idmnc(1));
% 
% fprintf(FFile, '%2.5f, ', FeaturesD.autoc(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.contr(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.corrm(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.corrp(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.cprom(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.cshad(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.dissi(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.energ(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.entro(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.homom(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.homop(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.maxpr(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.sosvh(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.savgh(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.svarh(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.senth(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.dvarh(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.denth(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.inf1h(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.inf2h(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.indnc(1));
% fprintf(FFile, '%2.5f, ', FeaturesD.idmnc(1));
% fprintf(FFile, '\n');
% fclose(FFile);
end





	