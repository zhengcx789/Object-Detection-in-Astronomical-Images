function BW_seedpoints=select_seedpoints(Global_BW,objects_num)
% This function is for selecting the large objects in the global detection results 
% as the seed points for the watershed method
%
% Input: im_ori = an astronomical image
%        Global_BW = a binary image which is objects detection results,the value 1 
%                    in it means that this pixel belongs to an object
%        objects_num: the number of selected seed points
% Output:
%       BW_seedpoints = a binary image(masks),the value 1 in it means that
%                       this pixel belongs to an seed points)  

BW2=imerode(Global_BW,strel('disk',8));
BW2=imdilate(BW2,strel('disk',8));
[BWlabeled,num2]=bwlabel(BW2,8);
BW_seedpoints=zeros(size(BW2));
for i=1:num2
    BWmask=(BWlabeled==i);
    subarea(i)=sum(BWmask(:));
end

[sortI,sortindx]= sort(subarea,'descend');
largeOBindx=sortindx(1:objects_num);
for i=1:size(largeOBindx,2)
    largeOBr=find(BWlabeled==largeOBindx(i));
    BW_seedpoints(largeOBr)=1;
end