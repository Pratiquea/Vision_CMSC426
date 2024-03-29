% MyRotobrush.m  - UMD CMSC426, Fall 2018
% This is the main script of your rotobrush project.
% We've included an outline of what you should be doing, and some helful visualizations.
% However many of the most important functions are left for you to implement.
% Feel free to modify this code as you see fit.

% Some parameters you need to tune:
WindowWidth = 40;  
ProbMaskThreshold = -1; 
NumWindows= 30; 
BoundaryWidth = 3;

% Load images:
fpath = '../frames/Frames1';
files = dir(fullfile(fpath, '*.jpg'));
imageNames = zeros(length(files),1);
images = cell(length(files),1);

for i=1:length(files)
    imageNames(i) = str2double(strtok(files(i).name,'.jpg'));
end

imageNames = sort(imageNames);
imageNames = num2str(imageNames);
imageNames = strcat(imageNames, '.jpg');

for i=1:length(files)
    images{i} = im2double(imread(fullfile(fpath, strip(imageNames(i,:)))));
end

% NOTE: to save time during development, you should save/load your mask rather than use ROIPoly every time.
mask = roipoly(images{1});
% Sss = load('turtle_mask.mat');
% mask = Sss.mask;

imshow(imoverlay(images{1}, boundarymask(mask,8),'red'));
% set(gca,'position',[0 0 1 1],'units','normalized')
% F = getframe(gcf);
% [I,~] = frame2im(F);
% imwrite(I, fullfi1le(fpath, strip(imageNames(1,:))));
% outputVideo = VideoWriter(fullfile(fpath,'video.mp4'),'MPEG-4');
% open(outputVideo);
% writeVideo(outputVideo,I);

% Sample local windows and initialize shape+color models:
[mask_outline, LocalWindows] = initLocalWindows(images{1},mask,NumWindows,WindowWidth,true);

ColorModels = ...
    initColorModels(images{1},mask,mask_outline,LocalWindows,BoundaryWidth,WindowWidth);

% You should set these parameters yourself:
fcutoff = 0.85;
SigmaMin = 2;
SigmaMax = WindowWidth;
R = 2;
A = (SigmaMax - SigmaMin)/((1 - fcutoff)^R);
ShapeConfidences = ...
    initShapeConfidences(mask_outline,LocalWindows,ColorModels,...
    WindowWidth, SigmaMin, A, fcutoff, R);

showColorConfidences(images{1}, mask_outline, ColorModels, LocalWindows, WindowWidth);

% Show initial local windows and output of the color model:
% imshow(images{1})
% hold on
% showLocalWindows(LocalWindows,WindowWidth,'r.');
% hold off
set(gca,'position',[0 0 1 1],'units','normalized')
F = getframe(gcf);
[I,~] = frame2im(F);

step_ = 5;
%%% MAIN LOOP %%%
% Process each frame in the video.
for prev=1:(length(files)-1)
    curr = prev+step_;
    fprintf('Current frame: %i\n', curr)
    
    %%% Global affine transform between previous and current frames:
    [warpedFrame, warpedMask, warpedMaskOutline, warpedLocalWindows] = calculateGlobalAffine(...
        images{prev}, images{curr}, mask, mask_outline,LocalWindows);
    %debug global affine
%     figure;
%     imshow(warpedFrame);
%     hold on
%     showLocalWindows(warpedLocalWindows,WindowWidth,'r.');
%     hold off
    %%% Calculate and apply local warping based on optical flow:
    NewLocalWindows = ...
        localFlowWarp(warpedFrame, images{curr}, warpedLocalWindows,warpedMask,WindowWidth);
    
    % Show windows before and after optical flow-based warp:
    imshow(images{curr});
    hold on
    showLocalWindows(warpedLocalWindows,WindowWidth,'r.');
    showLocalWindows(NewLocalWindows,WindowWidth,'b.');
    hold off
    
    %%% UPDATE SHAPE AND COLOR MODELS:
    % This is where most things happen.
    % Feel free to redefine this as several different functions if you prefer.
    [ ...
        mask, ...
        LocalWindows, ...
        ColorModels, ...
        ShapeConfidences, ...
    ] = ...
    updateModels(...
        NewLocalWindows, ...
        LocalWindows, ...
        images{curr}, ...
        warpedMask, ...
        warpedMaskOutline, ...
        WindowWidth, ...
        ColorModels, ...
        ShapeConfidences, ...
        ProbMaskThreshold, ...
        fcutoff, ...
        SigmaMin, ...
        R, ...
        A ...
    );

    mask_outline = bwperim(mask,4);

    % Write video frame:
    imshow(imoverlay(images{curr}, boundarymask(mask,8), 'red'));
    set(gca,'position',[0 0 1 1],'units','normalized')
    F = getframe(gcf);
    [I,~] = frame2im(F);
    imwrite(I, fullfile(fpath, strip(imageNames(curr,:))));
    writeVideo(outputVideo,I);

    imshow(images{curr})
    hold on
    showLocalWindXows(LocalWindows,WindowWidth,'r.');
    hold off
    set(gca,'position',[0 0 1 1],'units','normalized')
    F = getframe(gcf);
    [I,~] = frame2im(F);
end

close(outputVideo);
