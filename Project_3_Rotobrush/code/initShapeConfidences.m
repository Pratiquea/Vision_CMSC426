function ShapeConfidences = initShapeConfidences(MaskOutline, LocalWindows, ColorConfidences, WindowWidth, SigmaMin, A, fcutoff, R)
% INITSHAPECONFIDENCES Initialize shape confidences.  ShapeConfidences is a struct you should define yourself.

% fprintf('\nnumber of windows = ');
% length(LocalWindows)
% fprintf('\n');


% imshow(MaskOutline);
for win = 1:length(LocalWindows)
	% win
	% fprintf('\n');	
	fc = ColorConfidences(win).Confidences;
	if (fcutoff<fc && fc <=1);
		sigma_s = SigmaMin + A*(fc-fcutoff)^R;
	elseif(fc>=0 && fc<=fcutoff)
		sigma_s = SigmaMin;
	end
	% sigma_s
	x_c = LocalWindows(win,1);
	y_c = LocalWindows(win,2);
	window_maskoutline = MaskOutline( (y_c-WindowWidth/2):(y_c+WindowWidth/2),(x_c-WindowWidth/2):(x_c+WindowWidth/2) );
	% figure;
	% imshow(window_maskoutline);
	mask_elems = numel(window_maskoutline);

	% d = reshape(bwdist(MaskOutline( (y_c-WindowWidth/2):(y_c+WindowWidth/2),(x_c-WindowWidth/2):(x_c+WindowWidth/2) )),[mask_elems,1]);
	d = bwdist(window_maskoutline);
	ShapeConfidences(win).Confidences = 1-exp(-d.^2/sigma_s^2);
% 	figure(win);
% 	imshow(ShapeConfidences(win).Confidences);
end

end
