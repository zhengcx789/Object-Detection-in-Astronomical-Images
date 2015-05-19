function statistics=compt_statistics(im,objects_label)
% This function computes some statistics of detected objects, 
% such as center, area, magnitude and so on.
% Input:
%       im = the astronomical images
%       objects_label = the matrix labels the detected objects, different values larger 
%                       than 0 in this matrix represent different objects                    
% Output:
%       statistics = the statistics of detected objects

props=regionprops(objects_label,'Area'); % Regionprops is a function for computing the region properties
num=max(objects_label(:));

%Preallocating array
statistics(num) = charindex;

for x=1:size(im,1)
    for y=1:size(im,2)
    
        i=objects_label(x,y);
        
        if(objects_label(x,y)==0)
            continue;
        end
        if(any(statistics(i).area))
            continue;
        end
        statistics(i).area=props(i, 1).Area; % The areas of the detected objects

        %%%%%%% Estimate the range of the sweeping box window %%%%%%%%%%%%%
        max_w = round(sqrt(statistics(i).area)*3);
        max_h = round(sqrt(statistics(i).area)*3);

        min_bx = (x-max_w);
        min_by = (y-max_h);
        max_bx = (x+max_w);
        max_by = (y+max_h);

        min_bx=max([min_bx 1]);
        min_by=max([min_by 1]);
        max_bx=min([max_bx size(im,1)]);
        max_by=min([max_by size(im,2)]);

        o_label=objects_label(min_bx:max_bx,min_by:max_by);
        %imshow(o_label)
        
        %%%%%%%%%%%%%% Locate the pixels for the detected object (global px) %%%%%%%%%%%%
        z_i=(o_label==i);
        [r2,c2]=find(z_i==1);
        
        r=r2+min_bx-1; 
        c=c2+min_by-1;

        %%%%%%%%%%%%%% Compute the max/min X/Y cordinate of the detected object %%%%%%%%%%%%%%

        XMIN_IMAGE=min(c);% Minimum X-coordinate 
        YMIN_IMAGE=min(r);% Minimum Y-coordinate
        XMAX_IMAGE=max(c);% Maximum X-coordinate 
        YMAX_IMAGE=max(r);% Maximum Y-coordinate 
        statistics(i).MaxMinXY=[XMIN_IMAGE,XMAX_IMAGE,YMIN_IMAGE,YMAX_IMAGE];

        %fprintf('minmax done\n');

        %%%%%%%%%%%%%% Compute the mean intensity,magnitude of the detected object  %%%%%%%%%%%%%%
        si = im(YMIN_IMAGE:YMAX_IMAGE,XMIN_IMAGE:XMAX_IMAGE).*z_i(min(r2):max(r2),min(c2):max(c2));

        sumgray=sum(si(:)/1000);
        mgray=sumgray/size(r,1)*1000; % The mean intensity 
        statistics(i).meanIntensity=mgray;

        if sumgray>0;
            r_mag=0-2.5*log10(sumgray*1000); % The magnitude 
        else
            r_mag=99;  % Mark the magnitude which is a imaginary number as '99'.
        end
        statistics(i).r_mag=r_mag;

        %%%%%%%%%%%%%% Compute the barycenter and variance of the detected object  %%%%%%%%%%%%%%
        %sd=[];
        sumCx=0;
        sumCy=0;
        sumC2x=0;
        sumC2y=0;
        sumC2xy=0;

        % Factor in the shift of bounding box
        r2=r-YMIN_IMAGE+1; 
        c2=c-XMIN_IMAGE+1;

        for t=1:size(r2,1)
            %sd(t)=(si(r2(t),c2(t))-mgray)^2;
            sumCx=sumCx+si(r2(t),c2(t))/1000*c(t);
            sumCy=sumCy+si(r2(t),c2(t))/1000*r(t);
            sumC2x=sumC2x+si(r2(t),c2(t))/1000*c(t)*c(t);
            sumC2y=sumC2y+si(r2(t),c2(t))/1000*r(t)*r(t);
            sumC2xy=sumC2xy+si(r2(t),c2(t))/1000*c(t)*r(t);
        end

        X_IMAGE=sumCx*1000/sum(si(:));  % Weight centroid (X-axis)
        Y_IMAGE=sumCy*1000/sum(si(:));  % Weight centroid (Y-axis)
        X2_IMAGE=sumC2x*1000/sum(si(:))-X_IMAGE^2; % Variance along X-axis
        Y2_IMAGE=sumC2y*1000/sum(si(:))-Y_IMAGE^2; % Variance along Y-axis
        XY_IMAGE=sumC2xy*1000/sum(si(:))-X_IMAGE*Y_IMAGE;  % Covariance between X-axis and Y-axis
        statistics(i).barycenter=[X_IMAGE,Y_IMAGE];
        statistics(i).variance=[X2_IMAGE,Y2_IMAGE,XY_IMAGE];

        %%%%%%%%%%%%%% Compute the maximun and minimum intensities of the detected object %%%%%%%%%%%%%%
        mingrayvalue=inf;
        for i1=1:size(r,1)
            if  im(r(i1,1),c(i1,1))<mingrayvalue
                mingrayvalue=im(r(i1,1),c(i1,1));
            end
        end
        statistics(i).minIntensity= mingrayvalue; % The minimum intensity
        statistics(i).maxIntensity= max(si(:)); % The maximum intensity 

        %%%%%%%%%%%%%% Compute the major axis and minor axis of the detcted object%%%%%%%%%%%%%%
        re1=(X2_IMAGE+Y2_IMAGE)/2;
        re2=sqrt(((X2_IMAGE-Y2_IMAGE)/2)^2+XY_IMAGE^2);

        A_IMAGE=real(sqrt(re1+re2));   % The semi-major axis
        B_IMAGE=real(sqrt(re1-re2));   % The semi-minor axis

        if X2_IMAGE-Y2_IMAGE~=0
            tan2theta=2*XY_IMAGE/(X2_IMAGE-Y2_IMAGE);
            THETA_IMAGE=atan(tan2theta)/2;   % THETA_IMAGE is the position-angle between the major axis and the X image axis  (radians).
            if THETA_IMAGE>0 && XY_IMAGE<0
                THETA_IMAGE=THETA_IMAGE-pi/2;
            end
            if THETA_IMAGE<0 && XY_IMAGE>0
                THETA_IMAGE=THETA_IMAGE+pi/2;
            end
        else
            THETA_IMAGE=pi/4.0;
        end
        statistics(i).Smajoraxis=A_IMAGE;
        statistics(i).Sminoraxis=B_IMAGE;
        statistics(i).Theta= THETA_IMAGE;

    end
end
end


