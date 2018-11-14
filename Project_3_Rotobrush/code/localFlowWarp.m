function [NewLocalWindows] = localFlowWarp(WarpedPrevFrame, CurrentFrame, LocalWindows, Mask, Width)
% LOCALFLOWWARP Calculate local window movement based on optical flow between frames.
% TODO
opticFlow = opticalFlowFarneback;
% while win = 1:length(LocalWindows)
    
	I1 = im2double(rgb2gray(WarpedPrevFrame));
	I2 = im2double(rgb2gray(CurrentFrame));
	% opticalFlow = vision.OpticalFlow('ReferenceFrameSource', 'Input port')                                                          
	flow = estimateFlow(opticFlow, I1);%, I2);
	% opticalFlow = vision.OpticalFlow('ReferenceFrameSource', 'Input port', 'Method', 'Lucas-Kanade', 'OutputValue', 'Horizontal and Vertical Components in complex form')
 %    o=step(opticalFlow, I2, I1);

    % frameRGB = readFrame(vidReader);
    % frameGray = rgb2gray(frameRGB);

    % flow = estimateFlow(opticFlow,I1,I2); 

    imshow(CurrentFrame) 
    hold on
    plot(flow,'DecimationFactor',[5 5],'ScaleFactor',2)
    hold off 
% end

end

