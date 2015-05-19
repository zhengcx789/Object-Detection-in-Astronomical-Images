function Final_Label=label_objects(im,BW,FilterName,FilterSize,FilterSigmal,deblend_name,DEBLEND_MINCONT,DEBLEND_NTHRESH)
% This function to label each individual object by different value.
%
% Input: im = the astronomical image
%        BW = the binary image which is the objects detection results, and
%             the value 1 in it means that this pixel belongs to an object 
%        FilterName = the filter method
%        FilterSize = the size of filter window
%        FilterSigmal = the standard deviation SIGMA (positive)
%        deblend_name = the name of the deblending method
%        DEBLEND_MINCONT = the number of deblending subthresholds for the
%                          deblending function "deblend_Multi"
%        DEBLEND_NTHRESH = the minimum contrast parameter for deblending for the
%                          deblending function "deblend_Multi"
% Output:
%         Final_Label = the label matrix of all detected objects,different values 
%                       larger than 0 in this matrix represent different objects

[GlobalLabel,num]=bwlabel(BW,8);
OBJprops=regionprops(GlobalLabel,'Area');
Final_Label=GlobalLabel;
g = fspecial(FilterName,FilterSize,FilterSigmal);
img_gauss = conv2(im, g, 'same');
newlabel=num+1;
for nu=1:num
    if OBJprops(nu).Area>5
        [rind,cind]=find(GlobalLabel==nu);
        BWobjects=(GlobalLabel==nu);
        imlocal_BW=BWobjects(min(rind):max(rind),min(cind):max(cind));
        imlocal=img_gauss(min(rind):max(rind),min(cind):max(cind));
        
        % Select the deblend method
        if strcmp(deblend_name,'deblend_Multi')
            deblending_Label=deblend_Multi(imlocal,imlocal_BW,DEBLEND_MINCONT,DEBLEND_NTHRESH);
            
        elseif strcmp(deblend_name,'deblend_Wa')
            deblending_Label=deblend_Wa(imlocal,imlocal_BW);
            
        else
            disp('unknow deblend method')
        end
        
        numOutlabel=unique(deblending_Label);
        LabelNew=deblending_Label;
        if size(numOutlabel,1)>2
            sortnum=sort(numOutlabel,'ascend');
            LabelNew(deblending_Label==sortnum(2))=nu;
            for nuu=3:size(sortnum,1)
                Lnuu=sortnum(nuu);
                LabelNew(deblending_Label==Lnuu)=newlabel;
                newlabel=newlabel+1;
            end
            Final_Label(min(rind):max(rind),min(cind):max(cind))=Final_Label(min(rind):max(rind),min(cind):max(cind)).*(1-imlocal_BW)+LabelNew.*imlocal_BW;
        end
    end
end
end