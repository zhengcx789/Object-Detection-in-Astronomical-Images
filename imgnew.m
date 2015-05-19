function Imagenew=imgnew(img,Minvalue,Maxvalue)
% This function is for clipping the image "img" intensity to form "Imagenew"
% which only retains the intensity between the minimum value "Minvalue" and
% maximum value "Maxvalue" in the image "img".
%
% Input:
%       Img = Input image
%       Minvalue = the minimum value
%       Maxvalue = the maximun value
% Output:
%       Imagenew = An array which only retains the intensity between the Minvalue and Maxvalue

Imagenew=[];
p=1;
for i=1:size(img,1)
    for j=1:size(img,2)
        if img(i,j)>Minvalue && img(i,j)<Maxvalue
            Imagenew(p,1)=img(i,j);
            p=p+1;
        end
    end
end