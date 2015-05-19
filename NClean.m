function BWnew=NClean(im,Global_BW,FilterName,FilterSize,FilterSigmal,rincrease,Thresh)
% This function is for cleaning the artifacts(noise) around large objects
%
% Input: im = the astronomical image
%        Global_BW= the binary image of objects detection results in the
%                     image "im"(the value 1 in "imlocal_BW" represnts the objects,
%                     the value 0 represents the background)
%        FilterName = the filter method
%        FilterSize = the size of filter window
%        FilterSigmal = standard deviation SIGMA (positive)
%        rincrease = the parameter for expending the local area
%        Thresh = The parameters for defining how large the object need to 
%                 be processed (clean the noise around it)     
% Output:
%        BWnew = the results of removing the noise around large objects (A binary image) 


[labeled,num]=bwlabel(Global_BW,8);
characters=regionprops(labeled,'Area','BoundingBox'); % regionprops is the function of computing the region properties

%%%%%%%%%%%%%% Compute the mean area of detcted objects %%%%%%%%%%%%%%
pp=1;
area=[];
for i=1:size(characters,1)
    area(pp,1)=characters(i, 1).Area;
    pp=pp+1;
end
meanarea=mean(area);

%%%%%%%%%%%%%% Cleaning the artifacts(noise) around large objects %%%%%%%%%%%%%%
BWnew=Global_BW;
for nn=1:size(characters,1)
    if characters(nn,1).Area>Thresh*meanarea
        % First step: re-detection the area of large objects for noise suppression
        box=characters(nn,1).BoundingBox; 
        oldBWarea=Global_BW(max(ceil(box(1,2))- rincrease,1):min(floor(box(1,2)+box(1,4))+ rincrease,size(im,1)),max(ceil(box(1,1))- rincrease,1):min(floor(box(1,1)+box(1,3))+ rincrease,size(im,2)));
        sub_im=im(max(ceil(box(1,2))- rincrease,1):min(floor(box(1,2)+box(1,4))+ rincrease,size(im,1)),max(ceil(box(1,1))- rincrease,1):min(floor(box(1,1)+box(1,3))+ rincrease,size(im,2)));
        sub_normal=sqrt(mat2gray(sub_im));
        BG=BGsubtract(sub_normal);
        sub_BG=sub_normal-BG;
        sub_histeq=histeq(sub_BG);
        filter_g = fspecial(FilterName,FilterSize,FilterSigmal);
        sub_smooth = conv2(sub_histeq, filter_g, 'same');
        Thresh_BW = graythresh(sub_smooth);  
        newsubBW2= im2bw(sub_smooth,Thresh_BW);
        while sum(newsubBW2(:))-sum(oldBWarea(:))>sum(oldBWarea(:))/3
            Thresh_BW = Thresh_BW+0.01;
            newsubBW2= im2bw(sub_smooth,Thresh_BW);
        end
        
        % Second step: apply the open and close operator for cleaning noise around large objects
        newsubBW=imopen(newsubBW2,strel('disk',2));
        newsubBW=imclose(newsubBW,strel('disk',2));
        
        row_bottom= min(floor(box(1,2)+ box(1,4))+ rincrease,size(im,1))-2;
        row_up= max(ceil(box(1,2))- rincrease,1)+2;
        column_up=  max(ceil(box(1,1))- rincrease,1)+2;
        column_bottom=min(floor(box(1,1)+ box(1,3))+ rincrease,size(im,2))-2;
        BWnew(row_up:row_bottom,column_up:column_bottom)=newsubBW(3:size(newsubBW,1)-2,3:size(newsubBW,2)-2);
    end
end