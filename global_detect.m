function [Globaldetect_BW,Thresh_D]=global_detect(im,FilterName,FilterSize,FilterSigmal)
% This function detects all pixel set larger than a threshold in whole image as the objects 
% Input:
%      im = the astronomical image
%      FilterName = the filter method
%      FilterSize = the size of filter window
%      FilterSigmal = standard deviation SIGMA (positive)
% Output:
%      Globaldetect_BW = the detection results (a binary image) where the value 1 
%                        signifies a pixel belonging to an object
%      Thresh_D = the threshold for detecting objects

warning ('off');
%%%%%%%%%%%%%% Use filter to remove noise %%%%%%%%%%%%%%%%%
g = fspecial(FilterName,FilterSize,FilterSigmal);
img = conv2(im, g, 'same');

%%%%%%%%%%%%%% Gray scaling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imgnormal_m=mat2gray(img);
imgnormal_s=sqrt(imgnormal_m);

%%%%%%%%%%%%% Background subtraction %%%%%%%%%%%%%%%%%%%
BG=BGsubtract(imgnormal_s);
imgnormal=imgnormal_s-BG;

%%%%%%%%%%%%%% Histogram equalization %%%%%%%%%%%%%%%%%%%
histeqImg=histeq(imgnormal,256);

%%%%%%%%%%%%%% Image Binaryzation by Otsu's method %%%%%%%%%%%%%
Thresh_D = graythresh(histeqImg);  % Compute binarization threshold by using Otsu's method
Globaldetect_BW = im2bw(histeqImg,Thresh_D);

