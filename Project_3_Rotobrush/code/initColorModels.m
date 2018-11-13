function ColorModels = initializeColorModels(IMG, Mask, MaskOutline, LocalWindows, BoundaryWidth, WindowWidth)
% INITIALIZAECOLORMODELS Initialize color models.  ColorModels is a struct you should define yourself.
%
% Must define a field ColorModels.Confidences: a cell array of the color confidence map for each local window.



k = 1;
signma_c = WindowWidth/2;
IMG = rgb2lab(IMG);

%%%%%%%%%%%%%%%%%%%% GMM forground %%%%%%%%%%%%%%%%%%%%
data_f = [];
data_b = [];
for win = 1:length(LocalWindows)
	% fprintf('\n');
	x_c = LocalWindows(win,1);
	% fprintf(',');
	y_c = LocalWindows(win,2);
	% figure
	% imshow(Mask( (y_c-(WindowWidth/2)):(y_c+(WindowWidth/2)),(x_c-(WindowWidth/2)):(x_c+(WindowWidth/2)) ));
	%hold on
	%plot(xc,yc,'go');
	%hold off
	crop_img = IMG((y_c-(WindowWidth/2)):(y_c+(WindowWidth/2)),(x_c-(WindowWidth/2)):(x_c+(WindowWidth/2)),:);
	% figure
	% imshow(crop_img);
	[points_x_f,points_y_f] = find(Mask( (y_c-(WindowWidth/2)):(y_c+(WindowWidth/2)),(x_c-(WindowWidth/2)):(x_c+(WindowWidth/2)) ));
	[points_x_b,points_y_b] = find(~Mask( (y_c-(WindowWidth/2)):(y_c+(WindowWidth/2)),(x_c-(WindowWidth/2)):(x_c+(WindowWidth/2)) ));
	data_f = [];
	data_b = [];
	for channel = 1:3
		col_f = [];
		col_b = [];
		for pixel = 1:length(points_x_f)
			col_f = vertcat(col_f,crop_img(points_y_f(pixel),points_x_f(pixel),channel) );
		end
		data_f = [data_f col_f];

		for pixel = 1:length(points_x_b)
			col_b = vertcat(col_b, crop_img(points_y_b(pixel),points_x_b(pixel),channel) );
		end
		data_b = [data_b col_b];
	end
	% figure
	% imshow(IMG);
	% hold on
	% plot(points_x_b,points_y_b,'ro');
	% hold off
	% figure
	% plot3(data_f(:,1), data_f(:,2),data_f(:,3),'o');
	% figure
	% plot3(data_b(:,1), data_b(:,2),data_b(:,3),'o');
	gmm_f = fitgmdist(data_f,k,'RegularizationValue',0.01);
	gmm_b = fitgmdist(data_b,k,'RegularizationValue',0.01);
	elems1 = crop_img(:,:,1);
	elems2 = crop_img(:,:,2);
	elems3 = crop_img(:,:,3);
	length_elems = numel(elems1);
	col1 = reshape(elems1,[length_elems,1]);
	col2 = reshape(elems2,[length_elems,1]);
	col3 = reshape(elems3,[length_elems,1]);
	window_pixels = [col1 col2 col3];

	likelihood_f = pdf(gmm_f,window_pixels);
	likelihood_b = pdf(gmm_b,window_pixels);
	ColorModels(win).gmm_f = gmm_f;
	ColorModels(win).gmm_b = gmm_b;
	ColorModels(win).prob = likelihood_f./(likelihood_f+likelihood_b);
	mask_elems = numel(Mask( (y_c-WindowWidth/2):(y_c+WindowWidth/2),(x_c-WindowWidth/2):(x_c+WindowWidth/2) ));
	L_t = reshape(Mask( (y_c-WindowWidth/2):(y_c+WindowWidth/2),(x_c-WindowWidth/2):(x_c+WindowWidth/2) ),[mask_elems,1]);
	d = reshape(bwdist(MaskOutline( (y_c-WindowWidth/2):(y_c+WindowWidth/2),(x_c-WindowWidth/2):(x_c+WindowWidth/2) )),[mask_elems,1]);
	w_c = exp(-d.^2/signma_c^2); 
	ColorModels(win).Confidences =1-(sum( (L_t-ColorModels(win).prob).*w_c )/sum(w_c));
end


%each window
	%select foreground pixels
		%calculte gmm for foreground pixels
			%[OUTPUT = mu and sigma for each k]
	%select background pixels
		%calculte gmm for background pixels
			%[OUTPUT = mu and sigma for each k]
