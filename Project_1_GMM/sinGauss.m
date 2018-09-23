%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Using color thresholding %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%read all the files in the directory
path=dir('train_images/*.jpg');
nFiles=length(path);
%vert stack to store R,G,B channels of entire dataset
imageStack=zeros(1,3);
orange_pixels_lab = zeros(1,3);
imageStack_lab=zeros(1,3);
for i=1:nFiles-1
    %currImagePath=path(i).folder+path(i).name;
    currImagePath=fullfile(path(i).folder, path(i).name);
    %read the image
    currImage=imread(currImagePath);
    %converting image from RGB to l*a*b*   
    currImage_lab = rgb2lab(currImage);
    % Define thresholds for channel l for lab image
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

% visualize the orange pixel
% ll = 93;
% bb = length(orange_pixels_lab)/ll;
% visualize_orange_lab = reshape(orange_pixels_lab,ll,bb,3);
% imshow(visualize_orange_lab);
%converting the datatype to similar for all pixels
orange_pixels_lab = cast(orange_pixels_lab,'single');
% fprintf('dim of orange_pixels_lab = %d %d',size(orange_pixels_lab));


prior = 0.5;
cov_orange_lab = cov(orange_pixels_lab);
mean_orange_lab = mean(orange_pixels_lab);
threshold = 0.0000003;      %threhold for a* channel

distance_ball = zeros(1);
area_ball = zeros(1);
for images = 1:nFiles
    curr_path = fullfile(path(images).folder, path(images).name);
    [filepath,name,ext] = fileparts(curr_path);
    distance_ball = vertcat(distance_ball, str2num(name));

    i = imread(curr_path);
    i_lab = rgb2lab(i);
    % imshow(i_lab);
    maskedLABImage = i_lab;
    i_lab = im2double(i_lab);
    i_lab = reshape(i_lab,640*480,3);
    likelihood_orange_lab = mvnpdf([i_lab(:,1) i_lab(:,2) i_lab(:,3)], mean_orange_lab, cov_orange_lab);
    posterior_orange_lab = likelihood_orange_lab*prior;
    filtered_img = reshape(posterior_orange_lab,640,480);
    %setting all pixels that have a* and b* value greater than some threshold 'ta' and 'tb'
    sliderBW_lab = filtered_img>threshold;
      bw = sliderBW_lab;
    maskedLABImage(repmat(~sliderBW_lab,[1 1 3]))=0;
    % imshow(maskedLABImage);%,'InitialMagnification', 'fit');
    filename_i = strcat('masked_',name,'.jpg');
    imwrite(maskedLABImage,filename_i);
    
    stats = regionprops('table',bw,'Centroid',...
    'MajorAxisLength','MinorAxisLength');
    [val_1 ind_1] = max(stats.MajorAxisLength);
    [val_2 ind_2] = max(stats.MinorAxisLength);
    if (ind_1 ==ind_2)
        centers  = stats.Centroid(ind_1,:);
        radii = (stats.MajorAxisLength(ind_1)+stats.MinorAxisLength(ind_1))/4 ;
    else
        radii_1 = (stats.MajorAxisLength(ind_1)+stats.MinorAxisLength(ind_1))/4;
        radii_2 = (stats.MajorAxisLength(ind_2)+stats.MinorAxisLength(ind_2))/4;
        if radii_1>radii_2
            radii = radii_1;
            centers = stats.Centroid(ind_1,:);
        else
            radii = radii_2;
            centers = stats.Centroid(ind_2,:);
        end
    end
    area_ball = vertcat(area_ball, pi*(radii^2));
    imshow(sliderBW_lab)
    hold on
    viscircles(centers,radii);
    hold off  

    % pause(1);
end
distance_ball = distance_ball(2:length(distance_ball));
area_ball = area_ball(2:length(area_ball));

f=fit(area_ball, distance_ball,'poly4');
plot(f, area_ball, distance_ball);


%%%%%%%%%%%%%%%%parameters of model
p1 = 1.025e-10; 
p2 = -5.428e-07; 
p3 = 0.001034;
p4 = -0.8775;
p5 = 387.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% y = p1*x^4 + p2*x^3 + p3*x^2 + p4*x + p5;       

test_path = dir('test_images/*.jpg');
n = length(test_path);
for elem = 1:n
    curr_path = fullfile(test_path(elem).folder, test_path(elem).name);
    [filepath,name,ext] = fileparts(curr_path);
    % fprintf('distance ground truth for image %d is  = %d',elem,str2num(name));
    
    test_i = imread(curr_path);
    image_lab = rgb2lab(test_i);
    % imshow(i_lab);
    maskedLABImage = image_lab;
    image_lab = im2double(image_lab);
    image_lab = reshape(image_lab,640*480,3);
    likelihood_orange_lab = mvnpdf([image_lab(:,1) image_lab(:,2) image_lab(:,3)], mean_orange_lab, cov_orange_lab);
    posterior_orange_lab = likelihood_orange_lab*prior;
    filtered_img = reshape(posterior_orange_lab,640,480);
    %setting all pixels that have a* and b* value greater than some threshold 'ta' and 'tb'
    sliderBW_lab = filtered_img>threshold;
      bw = sliderBW_lab;
    maskedLABImage(repmat(~sliderBW_lab,[1 1 3]))=0;
    imshow(maskedLABImage);%,'InitialMagnification', 'fit');
    filename_i = strcat('masked_test_',name,'.jpg');
    imwrite(maskedLABImage,filename_i);
    % pause(1);
    stats = regionprops('table',bw,'Centroid',...
    'MajorAxisLength','MinorAxisLength');
    [val_1 ind_1] = max(stats.MajorAxisLength);
    [val_2 ind_2] = max(stats.MinorAxisLength);
    if (ind_1 ==ind_2)
        centers  = stats.Centroid(ind_1,:);
        radii = (stats.MajorAxisLength(ind_1)+stats.MinorAxisLength(ind_1))/4 ;
    else
        radii_1 = (stats.MajorAxisLength(ind_1)+stats.MinorAxisLength(ind_1))/4;
        radii_2 = (stats.MajorAxisLength(ind_2)+stats.MinorAxisLength(ind_2))/4;
        if radii_1>radii_2
            radii = radii_1;
            centers = stats.Centroid(ind_1,:);
        else
            radii = radii_2;
            centers = stats.Centroid(ind_2,:);
        end
    end
    test_area_ball = pi*(radii^2);
    imshow(sliderBW_lab)
    hold on
    viscircles(centers,radii);
    hold off  
    pause(1);
    test_distance_ball = p1*test_area_ball^4 + p2*test_area_ball^3 + p3*test_area_ball^2 + p4*test_area_ball + p5;
    fprintf('\ndistance for test image %d is found to be %d',elem, test_distance_ball);

    % pause(1);
end




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%code to test on single image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % curr_path = 'train_images/68.jpg';
    % i = imread(curr_path);
    %                 % fprintf('\nimage loaded');
    % i_lab = rgb2lab(i);
    % % imshow(i_lab);
    % maskedLABImage = i_lab;
    % i_lab = im2double(i_lab);
    % i_lab = reshape(i_lab,640*480,3);
    % likelihood_orange_lab = mvnpdf([i_lab(:,1) i_lab(:,2) i_lab(:,3)], mean_orange_lab, cov_orange_lab);
    % posterior_orange_lab = likelihood_orange_lab*prior;
    % filtered_img = reshape(posterior_orange_lab,640,480);
    % %setting all pixels that have a* and b* value greater than some threshold 'ta' and 'tb'
    % sliderBW_lab = filtered_img>threshold;
    % maskedLABImage(repmat(~sliderBW_lab,[1 1 3]))=0;
    % % imshow(maskedLABImage)
    % bw = sliderBW_lab;
    % bw2 = bwmorph(sliderBW_lab,'remove');
    % figure
    % imshow(bw2);\

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%code to plot circles by finding, center and radii
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % stats = regionprops('table',bw,'Centroid',...
    % 'MajorAxisLength','MinorAxisLength');
    % [val_1 ind_1] = max(stats.MajorAxisLength);
    % [val_2 ind_2] = max(stats.MinorAxisLength);
    % if (ind_1 ==ind_2)
    %     centers  = stats.Centroid(ind_1,:)
    %     radii = (stats.MajorAxisLength(ind_1)+stats.MinorAxisLength(ind_1))/4 
    % else
    %     radii_1 = (stats.MajorAxisLength(ind_1)+stats.MinorAxisLength(ind_1))/4;
    %     radii_2 = (stats.MajorAxisLength(ind_2)+stats.MinorAxisLength(ind_2))/4;
    %     if radii_1>radii_2
    %         radii = radii_1
    %         centers = stats.Centroid(ind_1,:)
    %     else
    %         radii = radii_2
    %         centers = stats.Centroid(ind_2,:)
    %     end
    % end
    % imshow(sliderBW_lab)
    % hold on
    % viscircles(centers,radii);
    % hold off    zeros(1)