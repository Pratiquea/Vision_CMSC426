%% functionname: function description
function my_harris(image_, window_size, corner_threshold)

	%defining some constants
	% window_size = 3;
	% corner_threshold = 7000000000;
	sigma = 1;
	k = 0.5;

	% image_ = imread('2.jpg');
	% figure(1);
	image_gray = rgb2gray(image_);
	% subplot(2,3,1);


	[dima, dimb] = size(image_gray);

	%derivaties operator in x and y direction
	dx = [-1 0 1;-1 0 1;-1 0 1];
	dy = [-1 -1 -1; 0 0 0; 1 1 1];

	% calculating x and y derivative
	Ix = conv2(double(image_gray),dx,'same');
	Iy = conv2(double(image_gray),dy,'same');
	% figure
	% imshow(Ix);
	% figure
	% imshow(Iy);
	% imwrite(Ix,'differential_in_x_direction.jpg')
	% imwrite(Iy,'differential_in_y_direction.jpg')

	%creating a gaussian filter window, which we will slide over the image
	gauss = fspecial('gaussian', window_size, sigma);

	%Calculating M matrix
	Ix_square = conv2(double(Ix.*Ix), gauss,'same');
	Iy_square = conv2(double(Iy.*Iy), gauss,'same');
	Ix_Iy = conv2(double(Ix.*Iy), gauss,'same');

	%Calculating R(Harris measure)
	%**** R = Det(M) - kTrace(M)^2 ****
	det_M = (Ix_square.*Iy_square - Ix_Iy.*Ix_Iy);
	Trace_M = k*((Ix_square+Iy_square).^2);
	R = det_M - Trace_M;
	R_r = (Ix_square.*Iy_square - Ix_Iy.*Ix_Iy)./(Ix_square + Iy_square + eps);
	R_invert = -R;
	max_R = max(max(R_invert));
	min_R = min(min(R_invert));
	ratio = 255/(max_R - min_R);
	R_remap = (R_invert - min_R)*ratio;

	R_thresh = R_invert.*(R_invert>corner_threshold);

	alpha = ((R_thresh))*0.7;

	% finding the local maxim for non maximum supression
	fx_order = 3;
	maxima = ordfilt2(R_r, fx_order^2,ones(fx_order));
	corners = (R_r == maxima) & (R_r >corner_threshold);
	[i,j] = find(corners);
	figure
	imshow(image_);
	hold on
	% title('my_harris corner detection');
	% op = (R_thresh > imdilate(R_thresh, [1 1 1; 1 0 1; 1 1 1]));
	[j,i] = find(R_thresh > imdilate(R_thresh, [1 1 1; 1 0 1; 1 1 1]));
	plot(i,j,'ys');
	title(sprintf('corner_threshold = %d', corner_threshold));
	% imshowpair(image_,op,'diff')
	% figure
	% imshow(output);X
end