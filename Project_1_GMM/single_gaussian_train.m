%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% ROI cropping using color thresholding %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%read all the files in the directory
path=dir('train_images/*.jpg');
nFiles=length(path);
mat_filename = 'mean_and_cov.mat';
%vert stack to store R,G,B channels of entire dataset
orange_pixels_lab = zeros(1,3);
imageStack_lab=zeros(1,3);
for i=1:nFiles-1
    %loading full path name 
    currImagePath=fullfile(path(i).folder, path(i).name);
    %read the image
    currImage=imread(currImagePath);
    %converting image from RGB to l*a*b*   
    currImage_lab = rgb2lab(currImage);
    
    % Define thresholds for channel l for lab image
    % these were determined using color threshold app in matlab
    channel_l_Min = 19.510;
    channel_l_Max = 91.394;
    % Define thresholds for channel a for lab image
    channel_a_Min = 14.923;
    channel_a_Max = 51.449;
    % Define thresholds for channel b for lab image
    channel_b_Min = 10.049;
    channel_b_Max = 46.717;

    % Create mask based on chosen histogram thresholds
    sliderBW_lab = (currImage_lab(:,:,1) >= channel_l_Min ) & (currImage_lab(:,:,1) <= channel_l_Max) & ...
                    (currImage_lab(:,:,2) >= channel_a_Min ) & (currImage_lab(:,:,2) <= channel_a_Max) & ...
                    (currImage_lab(:,:,3) >= channel_b_Min ) & (currImage_lab(:,:,3) <= channel_b_Max);
    BW_lab = sliderBW_lab;

    % Initialize output masked image based on input image.
    maskedRGBImage_lab = currImage_lab;
    % Set background pixels where BW is false to zero.
    maskedRGBImage_lab(repmat(~BW_lab,[1 1 3])) = 0;
    maskedRGBImage_lab=reshape(maskedRGBImage_lab,640*480,3);
    %gathering all orange pixels(lab) in an matrix
    for pixels = 1:(640*480)
        if BW_lab(pixels) == 1
            orange_pixels_lab = vertcat(orange_pixels_lab,maskedRGBImage_lab(pixels,:));
        end
    end
end    

orange_pixels_lab = cast(orange_pixels_lab,'single');

%calculating convergence value and mean from the above data set
cov_orange_lab = cov(orange_pixels_lab);
mean_orange_lab = mean(orange_pixels_lab);

save(mat_filename,'cov_orange_lab','mean_orange_lab');
fprintf('mat file saved');