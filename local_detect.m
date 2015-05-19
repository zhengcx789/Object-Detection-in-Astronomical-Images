function [Localdetect_BW,partition_W] =local_detect(im,seedPoints,sigmoid_slope,sigmoid_center_p,size_thresh,noise_T1,noise_T2,FilterName,FilterSize,FilterSigmal)
% This function divides the entire image into sub-regions and detects objects in each region.
% Input:
%       im = the input astronomical image
%       seedPoints = the banity image (matrix) that consists only of large objects as the
%                    seed points for watershed method (the value "1" represents the seed point)
%       sigmoid_slope = the slope of the sigmoid function
%       sigmoid_center_p = a parameter for computing the center of the sigmoid function: 
%                          center of sigmoid function = sigmoid_center_p * median(image intensity values)
%       size_thresh = the threshold for defining how large the object should
%                     be processed by the layer detection scheme
%       noise_T1,noise_T2 = the threshold for adjusting the filter window size,which should
%                           be adjusted based on the noise levels of all local regions of image.
%       FilterName = the filter method
%       FilterSize = the size of filter window
%       FilterSigmal = standard deviation SIGMA (positive)
% Output:
%       Localdetect_BW =  the result of local detection (the value 1 in the matrix
%                          represnts the objects, the value 0 represents the background)
%       partition_W = the sub-region division results,its different value represents different region

%%%%%%%%%%%%%% Irregular Sub-region Division %%%%%%%%%%%%%%%
warning ('off');
W_p=1-seedPoints; % Set the values in the position of seed points are 0,others are 1.
L=watershed(W_p); % Watershed transformation

partition_W=double(L);
partition_W=imdilate(partition_W, strel('disk',2));
L_zero=zeros(size(partition_W));
[rinx,cinx]=find(partition_W==0);
if ~isempty(rinx)
    L_zero(partition_W==0)=1;
end
[L_zero_label,labelnum]=bwlabel(L_zero,8);
for Li=1:labelnum
    lablearea=(L_zero_label==Li);
    [rind,cind]=find(lablearea==1);
    lablearea_dilate=imdilate(lablearea, strel('disk',1));
    [rind2,cind2]=find(lablearea_dilate==1);
    [LIc,LOCc] = ismember(cind2,cind);
    loc=find(LIc==0);
    newlabel=partition_W(rind2(loc(1)),cind2(loc(1)));
    L_zero(L_zero_label==Li)=newlabel;
end
partition_W=partition_W+L_zero; % partition_W is the division result

%%%%%%%%%%%%%%% Detect objects in each sub-region %%%%%%%%%%%%%%%
imghisteq=zeros(size(im));
imgsmooth=zeros(size(im));
Localdetect_BW=zeros(size(im));

for tt=1:max(partition_W(:))
    [indR2,indC2]=find(partition_W==tt);
    if ~isempty(indR2)
        lefttop=max(1,min(indR2)-2);
        leftbottom=min(max(indR2)+2,size(partition_W,1));
        righttop=max(1,min(indC2)-2);
        rightbottom=min(max(indC2)+2,size(partition_W,2));
        subregion=partition_W(lefttop:leftbottom,righttop:rightbottom);
        subregion(subregion~=tt)=0;
        subregion(subregion==tt)=1;
        mask_small=subregion;
        maskarea_size=(leftbottom-lefttop-2*2)*(rightbottom-righttop-2*2);
        img_local=im(lefttop:leftbottom,righttop:rightbottom);
        
        
        %%%%%%%%%%%%%% Grayscale stretching %%%%%%%%%%%%%%
        img_local_trans1=sqrt(mat2gray(img_local));
        sigmoid_center=min(sigmoid_center_p*median(img_local_trans1(:)),1);
        img_local_trans2= 1./(1+exp(-sigmoid_slope.*(double(img_local_trans1)-sigmoid_center)+eps));
        
        %%%%%%%%%%%%%% Background estimation %%%%%%%%%%%%%%%
        BG=BGsubtract(img_local_trans2);
        img_BGsub=img_local_trans2-BG;
        
        %%%%%%%%%%%%%% Histogram equalization %%%%%%%%%%%%%%%
        imghisteq_local=histeq(img_BGsub,256); % Histgram equalization on the local image
        imghisteq(lefttop:leftbottom,righttop:rightbottom)=imghisteq_local.*mask_small+imghisteq(lefttop:leftbottom,righttop:rightbottom);
        
        %%%%%%%%%%%%%%% Adaptive noise removal %%%%%%%%%%%%%%%
        noise_level(tt)= NoiseLevel(imghisteq_local); % Noise Level estimation
        if noise_level(tt)>noise_T1
            filterW=FilterSize;
        elseif noise_level(tt)<noise_T2
            filterW=FilterSize-4;
        else
            filterW=FilterSize-2;
        end
        
        g_filter = fspecial(FilterName,filterW, FilterSigmal);
        imgsmooth_local = conv2(imghisteq_local, g_filter, 'same');
        imgsmooth(lefttop:leftbottom,righttop:rightbottom)=imgsmooth_local.*mask_small+imgsmooth(lefttop:leftbottom,righttop:rightbottom);
        
        %%%%%%%%%%%%%%% Thresholding image by Otsu's method %%%%%%%%%%%%%%%
        level = graythresh(imgsmooth_local);
        imgBW_local = im2bw(imgsmooth_local ,level);
        
        %%%%%%%%%%%%%%% Using the layered detection to detect more faint objects %%%%%%%%%%%%%%%%%
        imgBW_local_old=imgBW_local;
        [BWlabel,num]=bwlabel(imgBW_local,8);
        props=regionprops(BWlabel,'Centroid','Area'); 
        Ob_area=[];
        for mn=1:size(props,1)
            Ob_area(mn)=props(mn,1).Area;
        end
        BWlabel=imdilate(BWlabel,strel('disk',5));
        img_BGsub2=img_BGsub;
        img_BGsub3=img_BGsub;
        
        areaind2=[];
        [areaind1,areaind2]=find(Ob_area>size_thresh); % Find the object area larger than "size_thresh" which should be processed by layer detection scheme
                
        if ~isempty(areaind2)
            for ari=1:size(areaind2,2)
                Ob_label(ari)=BWlabel(round(props(areaind2(ari),1).Centroid(1,2)),round(props(areaind2(ari),1).Centroid(1,1)));
                img_BGsub2(BWlabel==Ob_label(ari))=0;
            end
            
            for arii= 1:length(Ob_label)
                img_BGsub3(BWlabel==Ob_label(arii))=mean(img_BGsub2(:)); % Set the intesities of large objects to mean intensity of the background subtracted image without the large objects
            end
            
            imghisteq_local2=histeq(img_BGsub3,256); % Histgram equalization
            imghisteq(lefttop:leftbottom,righttop:rightbottom)=imghisteq_local2.*mask_small+imghisteq(lefttop:leftbottom,righttop:rightbottom);
            
            g_filter = fspecial(FilterName,filterW, FilterSigmal);
            imgsmooth_local2 = conv2(imghisteq_local2, g_filter, 'same'); % Smoothing image
            imgsmooth(lefttop:leftbottom,righttop:rightbottom)=imgsmooth_local2.*mask_small+imgsmooth(lefttop:leftbottom,righttop:rightbottom);
            
            Binary_Tresh = graythresh(imgsmooth_local2); % Compute the threshold by using Otsu's method
            imgBW_local2 = im2bw(imgsmooth_local2 ,Binary_Tresh);
            
            
            BWlarge=zeros(size(BWlabel));
            for ari=1:size(Ob_label,2)
                BWlarge(BWlabel==Ob_label(ari))=1; % The large objects whose size larger than "size_thresh"
            end
            
            %%%%%%%%%%%%%%%% Combine the detected faint objects and large objects %%%%%%%%%%%%%%%%%
            imgBW_local3=imgBW_local2+BWlarge;
            [BWlabel2,num]=bwlabel(imgBW_local3,8);
            BWlabel2_2=BWlabel2;
            
            for ari=1:size(areaind2,2)
                BWlabel2(BWlabel2==BWlabel2(round(props(areaind2(ari),1).Centroid(1,2)),round(props(areaind2(ari),1).Centroid(1,1))))=0; % The label matrix which only includes the small objects
            end
            
            BWlabel_new=zeros(size(BWlabel2));
            for ari=1:size(areaind2,2)
                BWlabel_new(BWlabel2_2==BWlabel2_2(round(props(areaind2(ari),1).Centroid(1,2)),round(props(areaind2(ari),1).Centroid(1,1))))=1; % The label matrix which only includes the large objects
            end
            
            BWlabelnew=imdilate(BWlabel_new,strel('disk',3));
            imgBW_local_C=BWlabelnew+BWlabel2; % Combine the large and small objects
            imgBW_local=(imgBW_local_C>0);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Deal with the extreme case where the binary threshold computed
        %%%% by Otsu's method is too small (this leads to detect too much noise)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        addArea_ratio(tt)=sum(sum(imgBW_local-imgBW_local_old))/(size(imgBW_local_old,1)*size(imgBW_local_old,2));
        addArea_ratio2(tt)=sum(sum(imgBW_local-imgBW_local_old))/sum(imgBW_local_old(:));
        
        th=0.01;
        whitepixelarea=sum(imgBW_local(:));
        if whitepixelarea > 0.8*maskarea_size || addArea_ratio(tt)>0.18 || addArea_ratio2(tt)>0.7
            Binary_Tresh = graythresh(imgsmooth_local);
            imgBW_local = im2bw(imgsmooth_local ,Binary_Tresh);
            whitepixelarea=sum(imgBW_local(:));
            while whitepixelarea > 0.8*maskarea_size
                Binary_Tresh = min(Binary_Tresh+th,1);
                imgBW_local = im2bw(imghisteq_local ,Binary_Tresh);
                th=th+0.01;
                whitepixelarea=sum(imgBW_local(:));
            end
        end
        Localdetect_BW(lefttop:leftbottom,righttop:rightbottom)=imgBW_local.*mask_small+Localdetect_BW(lefttop:leftbottom,righttop:rightbottom);
    end
end
end
