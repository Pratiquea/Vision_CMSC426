function [WarpedFrame, WarpedMask, WarpedMaskOutline, WarpedLocalWindows] = calculateGlobalAffine(IMG1,IMG2,Mask,MaskOutline,Windows)
% CALCULATEGLOBALAFFINE: finds affine transform between two frames, and applies it to frame1, the mask, and local windows.
IM1 = rgb2gray(IMG1);
IM2 = rgb2gray(IMG2);

%Debugging; show grayscale images
% figure
% subplot(1,2,1);
% imshow(IM1);
% subplot(1,2,2);
% imshow(IM2);

%detecting corner points
	%selecting only points that lie on object that we are tracking
	[mask_coordinates(:,1) mask_coordinates(:,2)] = find(Mask);
	min_x = min(mask_coordinates(:,1));
	max_x = max(mask_coordinates(:,1));
	min_y = min(mask_coordinates(:,2));
	max_y = max(mask_coordinates(:,2));
	allowance = 10;
	%debugging mask values
	% figure
	% imshow(IM1);
 	%rectangle('Position', [min_y-allowance min_x-allowance max_y-min_y+2*allowance max_x-min_x+2*allowance]);
 	ROI_rect = [min_y-allowance min_x-allowance max_y-min_y+2*allowance max_x-min_x+2*allowance];
	IM1_points = detectSURFFeatures(IM1,'MetricThreshold', 10, 'ROI',ROI_rect);
	IM2_points = detectSURFFeatures(IM2,'MetricThreshold', 10, 'ROI',ROI_rect);
	
	%Extract the features.
	[f1,vpts1] = extractFeatures(IM1,IM1_points);
	[f2,vpts2] = extractFeatures(IM2,IM2_points);

	%Retrieve the locations of matched points.
	indexPairs = matchFeatures(f1,f2) ;
	matchedPoints1 = vpts1(indexPairs(:,1));
	matchedPoints2 = vpts2(indexPairs(:,2));
	%debug surf features
	% figure; showMatchedFeatures(IM1,IM2,matchedPoints1,matchedPoints2);
	% legend('matched points 1','matched points 2');

	%estimating affine transform
	[tform,inlierPtsDistorted,inlierPtsOriginal] = ...
    estimateGeometricTransform(matchedPoints1,matchedPoints2,'affine');

	% Recover the original image from the distorted image.
	outputView = imref2d(size(IM1));
	IM1_warped = imwarp(IM2,tform,'OutputView',outputView);
	outputView = imref2d(size(IMG1));
	IMG1_warped = imwarp(IMG2,tform,'OutputView',outputView);
	WarpedFrame = IMG1_warped;
	WarpedMask = imwarp(Mask,tform,'OutputView',outputView);
	figure;
	imshowpair(WarpedMask,Mask);
	WarpedMaskOutline = bwperim(WarpedMask,4);
	%debugging; original and warped image montage
	figure;
	imshowpair(WarpedMaskOutline,MaskOutline);
	% figure; 
	% imshowpair(IMG1,IMG1_warped,'montage'); 
	% title('Warped image');

	%

%

end

