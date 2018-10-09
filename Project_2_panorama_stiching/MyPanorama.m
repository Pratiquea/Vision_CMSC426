path_ = dir('images/input/*.jpg');
num_images = length(path_);
image_tensor = [];
for image_ = 1:num_images
	currImagePath=fullfile(path_(image_).folder, path_(image_).name);
	current_image = rgb2gray(imread(currImagePath));
	image_tensor = cat(3,image_tensor, current_image);

	% subplot(1,2,1);
	% subplot(1,2,2);
	corner_score= cornermetric(image_tensor(:,:,image_));
	% subplot(1,2,1);
	% imagesc(corner_score);
	% title('cornermetric');
	local_maxima = imregionalmax(corner_score);

	N_best = 200;
	N_strong= sum(local_maxima(:)==1);
	[y x] = find(local_maxima);

	figure;
	imshow(image_tensor(:,:,image_));
	hold on
	plot(x,y,'rs');
	title('Output after imregionalmax');
	hold off
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%% Adaptive Non-Maxima Supression step %%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	anms_output = [];
	Euclidean_Distance = 0;
	for i = 1:N_strong
		anms_output(i).radial_dist = inf;
		anms_output(i).x = 0;
		anms_output(i).y = 0;
	end
	% Euclidean_Distance = 0;
	for i = 1:N_strong
		for j = 1:N_strong
			if (corner_score(y(j),x(j)) > corner_score(y(i),x(i)))
				Euclidean_Distance = cal_Euclidean_Distance(x(i),y(i),x(j),y(j));
			end
			if(Euclidean_Distance < anms_output(i).radial_dist)
				anms_output(i).radial_dist = Euclidean_Distance;
				anms_output(i).x = x(j);
				anms_output(i).y = y(j);
			end
		end
	end

	%%%%%%%%%%%%% Sorting array %%%%%%%%%%%%%%%
		%converting structure to cell
		cells = struct2cell(anms_output);
		%selecting the field to sort
		sortvals = cells(1,1,:);
		% converting cell to matrix
		mat = (cell2mat(sortvals));
		% sorting matrix
		[sorted_,ix] = sort(mat,'descend');
		% sorting the original structure according to the field(or matrix)
		anms_output_sorted = anms_output(ix);

		clear x1;
		clear y1;
		x1 = zeros(N_best,num_images);
		y1 = zeros(N_best,num_images);

	% selecting first N_best elements from the list(ones having the largest distance)
	for i = 1:N_best
		x1(i,image_) = anms_output_sorted(i).x;
		y1(i,image_) = anms_output_sorted(i).y;
	end

	figure;
	imshow(image_tensor(:,:,image_));
	hold on
	plot(x1(:,image_),y1(:,image_),'rs');
	title('After ANMS');
	hold off

	pad = 20;
	patch_size = 40;
	sigma_ = 1.2;
	features = zeros(((patch_size/4)+1)^2,length(x1(:,image_)),num_images);
	x_padded(:,image_) = x1(:,image_)+(patch_size/2);
	y_padded(:,image_) = y1(:,image_)+(patch_size/2);
	padded_image = padarray(image_tensor(:,:,image_),[pad pad],'both','symmetric');
	for i = 1:length(x1(:,image_))
		% fprintf('\nx1-x2 = %d-%d',x_padded(i)-(patch_size/2),x_padded(i)+(patch_size/2));
		patch_ = padded_image(y_padded(i)-(patch_size/2):y_padded(i)+(patch_size/2),x_padded(i)-(patch_size/2):x_padded(i)+(patch_size/2));
		filtered_patch = imgaussfilt(patch_,sigma_);
		downsize_patch = imresize(filtered_patch,0.25,'nearest');
		features(:,i,image_) = imresize(downsize_patch,[numel(downsize_patch),1]);
	end

	mu = mean(features(:,:,image_));
	sd = std(features(:,:,image_));
	norm_features(:,:,image_) = (features(:,:,image_)-mu)./sd;

	% plot last patch
	% figure
	% subplot(2,1,1)
	% imshow(patch_);
	% subplot(2,1,2);
	% imshow(filtered_patch);
end

best_match = ones(length(norm_features(1,:,1)),1,num_images-1)*inf;
second_best_match =  ones(length(norm_features(1,:,1)),1,num_images-1)*inf;
for image_ = 1:(num_images-1)
	for i = 1:length(norm_features(1,:,1))
		for j = 1:length(norm_features(1,:,1))
			X = norm_features(:,i,image_)-norm_features(:,j,image_+1);
			ssd = sum(X(:).^2);
			if(best_match(i,1,image_)>ssd)
				second_best_match(i,1,image_) = best_match(i,1,image_);
				best_match(i,1,image_) = ssd;
			elseif (second_best_match(i,1,image_)>ssd)
				second_best_match(i,1,image_) = ssd;		
			end
			Intensity_ssd(i,j,image_)  = ssd;
		end
	end
end


ssd_ratio = best_match./second_best_match;
prefered_pts = ssd_ratio>0.8;


%gauss = fspecial('gaussian', window_size, sigma);

function dist = cal_Euclidean_Distance(x,y,a,b)
	dist = (x-a)^2+(y-b)^2;
end
