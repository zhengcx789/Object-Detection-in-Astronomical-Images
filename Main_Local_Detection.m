% This is a script to execute the local objects detection method
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
FilterSigmal=1.5;                       % Standard deviation SIGMA (positive)
DEBLEND_MINCONT=0.001;                  % Number of deblending subthresholds
DEBLEND_NTHRESH=32;                     % Minimum contrast parameter for deblending
sigmoid_slope=40;                       % Slope of the sigmoid function
sigmoid_center_p=1;                     % Parameter for computing the center of the sigmoid function
objects_num=30;                         % Number of selected seed points (watershed)
Noise_T1=0.08;                          % Threshold for adjusting the filter window size
Noise_T2=0.05;                          % Threshold for adjusting the filter window size
Size_thresh=3000;                       % Threshold defining object size needed to be 
                                        %    process by layered detection scheme
Deblend_method='deblend_Multi';         % Deblending type
ShowImage=1;                            % Display figures during computation
Verbose=1;                              % Enable console output during computation


% Set the path for saving results
path_name =strcat('Local method detection results\');
if ~exist(path_name)
    mkdir(path_name);
end
path = [pwd,strcat('\',path_name,'\')];

% Read an astronomical image
ori_im=fitsread('data\R.fits');
if ShowImage figure,imshow(ori_im,[]),title('original image'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   Following is the LOCAL detection pipeline  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Detect the objects in entire image %%%%%%%%%%%%%%
if Verbose fprintf('Beginning local object detection...'); end
[G_BW,Thresh_G]=global_detect(ori_im,FilterName,FilterSize,FilterSigmal);
if Verbose fprintf('..Done!\n'); end

%%%%%%%%%%%%%% Select seed points for the watershed method %%%%%%%%%%%%%%
if Verbose fprintf('Selecting seed points for the watershed method...'); end
BW_seedpoints=select_seedpoints(G_BW,objects_num);
if Verbose fprintf('..Done!\n'); end

%%%%%%%%%%%%%% Local detection based on watershed region %%%%%%%%%%%%%%
if Verbose fprintf('Begining local detection method...'); end
[localdetect_BW,Subregion_label_W]=local_detect(ori_im,BW_seedpoints,sigmoid_slope,sigmoid_center_p,Size_thresh,Noise_T1,Noise_T2,FilterName,FilterSize,FilterSigmal);
if Verbose fprintf('..Done!\n'); end

%Show and save sub-region partition by watershed method
Lrgb2 = label2rgb(Subregion_label_W, 'jet');
figure,imshow(BW_seedpoints);hold on; % show the seedpoint in each region
himage = imshow(Lrgb2);
set(himage, 'AlphaData', 0.3);
title('sub-region division by Watershed')
% print(gcf,'-dbmp',[path,'sub_region division by Watershed.png']);


%%%%%%%%%%%%%% Use the median filter to remove noise %%%%%%%%%%%%%%
if Verbose fprintf('Removing noise with a median filter...'); end
Final_detection_result=medfilt2(localdetect_BW,[3,3]); 
if Verbose fprintf('..Done!\n'); end
% Show and save the final detection results 
if ShowImage figure, imshow(Final_detection_result,[]),title('Final detection result');
save([path,'Final_detection_result.mat'],'Final_detection_result');
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
if Verbose fprintf('..Done!\n'); end
% Save Final_statistics
save([path,'Final_statistics.mat'],'Final_statistics');
