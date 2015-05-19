% This is a script to execute the global objects detection method
% The output detection results "Final_detection_result" is a binary image,
% where the white pixels represent the detected objects.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors:Caixia Zheng, Jesus Pulido, Paul Thorman, Bernd Hamann
% Copyright(c) 2015, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

%%%%%%%%%%%%%% Set the parameters %%%%%%%%%%%%%%
FilterName='gaussian';                  % 2D filter type 
FilterSize=[7,7];                       % Kernel filter size
FilterSigmal=1.5;                       % standard deviation SIGMA (positive)
Rincrease=30;                           % Parameter controlling local area expansion
Size_scalor=2;                          % Minimum pixel object size for detection
DEBLEND_MINCONT=0.001;                  % Number of deblending subthresholds
DEBLEND_NTHRESH=32;                     % Minimum contrast parameter for deblending
Deblend_method='deblend_Wa';            % Deblending type
ShowImage=1;                            % Display figures during computation
Verbose=1;                              % Enable console output during computation

% Set the path for saving results
path_name =strcat('global method detection results\'); 
if ~exist(path_name)
    mkdir(path_name);
end
path = [pwd,strcat('\',path_name,'\')];  

% Read an astronomical image
ori_im=fitsread('data\Deep_32.fits'); 
if ShowImage figure,imshow(ori_im,[]),title('original image'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   Following is the GLOBAL detection pipeline  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Detect the objects in entire image %%%%%%%%%%%%%%
if Verbose fprintf('Beginning global object detection...'); end
[G_BW,BW_Thresh]=global_detect(ori_im,FilterName,FilterSize,FilterSigmal);
if Verbose fprintf('..Done!\n'); end

%%%%%%%%%%%%%% Clean the artifacts(noise) around large objects %%%%%%%%%%%%%%
if Verbose fprintf('Beginning artifact cleaning...'); end
Final_detection_result=NClean(ori_im,G_BW,FilterName,FilterSize,FilterSigmal,Rincrease,Size_scalor);
if Verbose fprintf('..Done!\n'); end
% Show  and save the final detection results
save([path,'Final_detection_result.mat'],'Final_detection_result');
if ShowImage figure,imshow(Final_detection_result); title('Final detection result');
print(gcf,'-dpng',[path,'Final_detection_result.png']);
end
 
%%%%%%%%%%%%%%% Deblend and label the detected objects %%%%%%%%%%%%%%%
if Verbose fprintf('Beginning deblending...'); end
Final_Label=label_objects(ori_im,Final_detection_result,FilterName,FilterSize,FilterSigmal,Deblend_method,DEBLEND_MINCONT,DEBLEND_NTHRESH);
if Verbose fprintf('..Done!\n'); end
% Save Final_Label
save([path,'Final_Label.mat'],'Final_Label')

%%%%%%%%%%%%%%% Compute some statistics of detected objects %%%%%%%%%%%%%%%
if Verbose fprintf('Computing final statistics...'); end
Final_statistics=compt_statistics(ori_im,Final_Label);
if Verbose fprintf('..Done!\n');end
% Save Final_statistics
save([path,'Final_statistics.mat'],'Final_statistics');
