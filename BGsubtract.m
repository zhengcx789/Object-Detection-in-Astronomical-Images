function BG=BGsubtract(img)
% This function estimates the background of an astronomical image
% Input:
%       img = the image
% Output:
%       BG = the estimated background 

sigmal_ori=std(img(:));
minvalue=median(img(:))-3*sigmal_ori;
maxvalue=median(img(:))+3*sigmal_ori;

newim=imgnew(img,minvalue,maxvalue); % Discard the image value which is not within median ¡Ànum*sigmal_ori

sigmal_new=std(newim);
median_new=median(newim);
minvalue_new=median_new-3*sigmal_new;
maxvalue_new=median_new+3*sigmal_new;
sigmal_diff=[];
sigmal_diff(1)=(sigmal_ori-sigmal_new)/sigmal_ori;
m=2;

if sigmal_ori==0
    BG=zeros(size(img,1),size(img,2));
else
    while max(newim<=minvalue_new)==1 || max(newim>=maxvalue_new)==1   % Repeat until all the remaining pixel values are within median ¡Ànum*¦Ò
        sigmal=sigmal_new;
        newim=imgnew(newim,minvalue_new,maxvalue_new);
        sigmal_new=std(newim);
        median_new=median(newim);
        minvalue_new=median_new-3*sigmal_new;
        maxvalue_new=median_new+3*sigmal_new;
        sigmal_diff(m)=(sigmal-sigmal_new)/sigmal;
        m=m+1;
    end
end

if max(sigmal_diff)<0.2
    BG=mean(newim);
else
    BG=2.5*median(newim)-1.5*mean(newim);
end
