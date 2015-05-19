function out_label2=deblend_Wa(imlocal,imlocal_BW)
% This function is for deblending the touching or overlapping objects
% directly using the watershed method
%
% Input: imlocal = the local smoothed image (it is a local area including 
%                  a object which needs to be checked and applied the deblending processing)
%        imlocal_BW = the local binary image(local masks) correspoinding to the same region of "imlocal",
%                     the value 1 in it means that this pixel belongs to an object)                      
% Output:
%          outLabel2 = the local label matrix of deblended objects,different values 
%                      larger than 0 in this matrix represent different objects  


out_label2=imlocal_BW;

% Find the boundary points of the detcted object
[start_i,start_j]=find(imlocal_BW==1);
characters_new.boundary = bwtraceboundary(imlocal_BW,[start_i(1),start_j(1)],'w',8,inf,'clockwise');
contour=characters_new.boundary;
contour_points=contour(1:size(contour,1)-1,:);

tag=0;
% Check the object weather can be fitted by an ellipse or not
if size(contour_points,1)>2 && size(unique(contour_points(:,1)),1)>1 && size(unique(contour_points(:,2)),1)>1  
    re_para(1) = (contour_points(1,2)-contour_points(2,2))/(contour_points(1,1)-contour_points(2,1));
    re_para(2) = contour_points(1,2) - re_para(1) *contour_points(1,1);
    for tt=3:size(contour_points,1)
        if contour_points(tt,2)~=re_para(1)*contour_points(tt,1)+re_para(2); 
            tag=1;
            break;
        end
    end
    
    
    if tag==1
        % Fit ellipse to the detected object 
        params = fitellipse(contour_points(:,2),contour_points(:,1)); %fit ellipse
        points_p=[contour_points(:,2),contour_points(:,1)];
        ellipse_shape_param=[params(1),params(2),max(params(3:4)),min(params(3:4)),params(5)];
        
        % Compute the distance between the boundary points of the detcted object and the fitting ellipse
        [Dist,RSS,XYproj,di] = Residuals_ellipse(points_p,ellipse_shape_param);       
        Distnorm2=norm(Dist,2);
        
        if Distnorm2>12  % Check the object weather needs to be deblended or not
            % Deblend the object using the watershed method
            mask_em = imextendedmax(imlocal,30);
            bw=imlocal_BW;
            D=-bwdist(~bw); 
            D=imimposemin(D,mask_em); 
            L1=watershed(D);
            out_label2=double(L1).*double(bw);
        end
    end
end