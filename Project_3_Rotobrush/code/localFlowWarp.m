function [NewLocalWindows] = localFlowWarp(WarpedPrevFrame, CurrentFrame, LocalWindows, Mask, Width)
% LOCALFLOWWARP Calculate local window movement based on optical flow between frames.
% TODO
NewLocalWindows = [];
opticFlow = opticalFlowFarneback;
   
	warped_prev_frame = im2double(rgb2gray(WarpedPrevFrame));
	current_frame = im2double(rgb2gray(CurrentFrame));
	% opticalFlow = vision.OpticalFlow('ReferenceFrameSource', 'Input port')                                                          
	flow = estimateFlow(opticFlow, warped_prev_frame);
    flow = estimateFlow(opticFlow, current_frame);
    [m n] = size(warped_prev_frame);
    %debugging: optical flow comapare
      figure;
      imshowpair(WarpedPrevFrame, CurrentFrame) ;
      hold on;
      plot(flow,'DecimationFactor',[5 5],'ScaleFactor',2)
      hold off ;
    for win = 1:length(LocalWindows)
        x_c = LocalWindows(win,1);
		y_c = LocalWindows(win,2);
        x_win = round([(x_c-(Width/2)):(x_c+(Width/2))]);
        y_win = round([(y_c-(Width/2)):(y_c+(Width/2))]);

        window_ = warped_prev_frame;
        for q = 1:m
            for p =1:n
                if( find(x_win == p) & find(y_win == q) )
                    window_(q,p) = 1;
                else
                    window_(q,p) = 0;
                end
            end
        end
%         imshow(window_);
        mask_2 = window_.*Mask;
%         flow2.Vx = flow.Vx.*mask_2;
%         flow2.Vy = flow.Vy.*mask_2;
%         flow2.Orientation = flow.Orientation.*mask_2;
%         flow2.Magnitude = flow.Magnitude.*mask_2;
        [x_flow, y_flow] = find(mask_2);
        v_x = 0;
        v_y = 0;
        len_ = length(y_flow);
        for q = 1:len_
            v_x = v_x + flow.Vx(x_flow(q),y_flow(q));
            v_y = v_y + flow.Vy(x_flow(q),y_flow(q));
        end
        v_x = v_x/len_;
        v_y = v_y/len_;
        x_star = x_c + v_x;
        y_star = y_c + v_y;
        NewLocalWindows(win,1) = x_star;
        NewLocalWindows(win,2) = y_star;
%         flow_total = flow_total+flow2;
    end
    %degugging; new flow vector visualize
    min_x = min(x_flow);
    max_x = max(x_flow);
    min_y = min(y_flow);
    max_y = max(y_flow);
    figure;
    imshow(warped_prev_frame);
    hold on
  	ROI_rect = [min_y min_x max_y-min_y max_x-min_x];
    rectangle('Position', ROI_rect,'EdgeColor','b');
    ROI_rect = [(x_c-(Width/2)) (y_c-(Width/2)) Width Width];
    rectangle('Position', ROI_rect,'EdgeColor','r');
    hold off
%     
%     figure;
%     imshow(current_frame);
%     hold on
%     showLocalWindows(NewLocalWindows, Width,'r.');
%     showLocalWindows(LocalWindows, Width,'b.');
%     hold off
end

