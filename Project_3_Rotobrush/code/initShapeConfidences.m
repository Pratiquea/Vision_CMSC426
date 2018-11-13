function ShapeConfidences = initShapeConfidences(MaskOutline, LocalWindows, ColorConfidences, WindowWidth, SigmaMin, A, fcutoff, R)
% INITSHAPECONFIDENCES Initialize shape confidences.  ShapeConfidences is a struct you should define yourself.

fprintf('\nnumber of windows = ');
length(LocalWindows)
fprintf('\n');

for win = 1:length(LocalWindows)
	win
	fprintf('\n');	
	fc = ColorConfidences(win).Confidences;
	if (fcutoff<fc && fc <=1);
		sigma_s = SigmaMin + A*(fc-fcutoff)^R;
	elseif(fc>=0 && fc<=fcutoff)
		sigma_s = SigmaMin;
	
	x_c = LocalWindows(win,1);
	y_c = LocalWindows(win,2);
	
	mask_elems = numel(MaskOutline( (y_c-WindowWidth/2):(y_c+WindowWidth/2),(x_c-WindowWidth/2):(x_c+WindowWidth/2) ));

	% d = reshape(bwdist(MaskOutline( (y_c-WindowWidth/2):(y_c+WindowWidth/2),(x_c-WindowWidth/2):(x_c+WindowWidth/2) )),[mask_elems,1]);
	d = bwdist(MaskOutline( (y_c-WindowWidth/2):(y_c+WindowWidth/2),(x_c-WindowWidth/2):(x_c+WindowWidth/2) ));
	ShapeConfidences(win).Confidences = 1-exp(-d.^2/sigma_s^2);
	figure
	imshow(ShapeConfidences(win).Confidences);
end

end
