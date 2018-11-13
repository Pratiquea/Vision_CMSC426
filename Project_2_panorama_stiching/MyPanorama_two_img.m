clear;
path_ = dir('images/input/Set1/*.jpg');
% images/input/Set1
% num_images = length(path_);
image_tensor = [];
image_1 = rgb2gray(imread('1.jpg'));
image_2 = rgb2gray(imread('2.jpg'));
image_tensor = cat(3,image_tensor, image_1);
image_tensor = cat(3,image_tensor, image_2);

num_images = 2;
N_best = 50;
ssd_ratio_threshold = 0.6;
best_match_val_threshold = 20;
x1 = zeros(N_best,num_images);
y1 = zeros(N_best,num_images);
for image_ = 1:num_images
	% currImagePath=fullfile(path_(image_).folder, path_(image_).name);
	% current_image = rgb2gray(imread(currImagePath));
	% image_tensor = cat(3,image_tensor, current_image);

	% subplot(1,2,1);
	% subplot(1,2,2);
	corner_score= cornermetric(image_tensor(:,:,image_));
	% subplot(1,2,1);
	% imagesc(corner_score);
	% title('cornermetric');
	local_maxima = imregionalmax(corner_score);


	N_strong= sum(local_maxima(:)==1);
	[y x] = find(local_maxima);

	% figure;
	% % subplot(1,2,1);
	% title('Output after imregionalmax');
	% imshow(image_tensor(:,:,image_));
	% hold on
	% plot(x,y,'rs');
	% hold off
	% pause(1);
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



	% selecting first N_best elements from the list(ones having the largest distance)
	for i = 1:N_best
		x1(i,image_) = anms_output_sorted(i).x;
		y1(i,image_) = anms_output_sorted(i).y;
	end
	% pause(1);
	% figure
	% subplot(1,2,2);
	% title('After ANMS');
	% imshow(image_tensor(:,:,image_));
	% hold on
	% pause(1);
	% plot(x1(:,image_),y1(:,image_),'rs');
	% hold off

	pad = 20;
	patch_size = 40;
	sigma_ = 2;
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

% best_match_val = ones(length(norm_features(1,:,1)),1,num_images-1)*inf;
% best_match_index = ones(length(norm_features(1,:,1)),1,num_images-1)*NaN;
% second_best_match_val =  ones(length(norm_features(1,:,1)),1,num_images-1)*inf;

% for image_ = 1:(num_images-1)
% 	for i = 1:length(norm_features(1,:,1))
% 		for j = 1:length(norm_features(1,:,1))
% 			X = norm_features(:,i,image_)-norm_features(:,j,image_+1);
% 			ssd = sum(X(:).^2);
% 			if(best_match_val(i,1,image_)>ssd)
% 				second_best_match_val(i,1,image_) = best_match_val(i,1,image_);
% 				second_best_match_index(i,1,image_)= best_match_index(i,1,image_);

% 				best_match_val(i,1,image_) = ssd;
% 				best_match_index(i,1,image_) = j;
% 			elseif (second_best_match_val(i,1,image_)>ssd)
% 				second_best_match_val(i,1,image_) = ssd;
% 				second_best_match_index(i,1,image_) = j;
% 			end
% 			Intensity_ssd(i,j,image_)  = ssd;
% 		end
% 	end
% end


% ssd_ratio = best_match_val./second_best_match_val;
% prefered_pts = (ssd_ratio<ssd_ratio_threshold)&(best_match_val<best_match_val_threshold);

% %sum(prefered_pts(:,1));
% best_match_index(~prefered_pts) = NaN;
% len_vec = [];
% for image_ = 1:num_images-1
% 	len_vec = vertcat(len_vec,sum(prefered_pts(:,image_)));
% end
% max_vec_len = max(len_vec);
% match_pt1_tensor = [];			%NaN*ones(max_vec_len,2,num_images-1);
% match_pt2_tensor = [];		%NaN*ones(max_vec_len,2,num_images-1);
% for image_ = 1:num_images-1
% 	match_pt1 = [];
% 	match_pt2 = [];
% 	% image_
% 	for i = 1:length(best_match_index(:,1,1))
% 		if(~isnan(best_match_index(i,1,image_)))
% 			match_pt1 = vertcat(match_pt1,[x1(i,image_),y1(i,image_)]);
% 			match_pt2 = vertcat(match_pt2,[x1(best_match_index(i,1,image_),image_+1),y1(best_match_index(i,1,image_),image_+1)]);
% 		end		
% 	end
% 	% length(match_pt1)
% 	if(length(match_pt1)<max_vec_len)
% 		for m = length(match_pt1)+1 : max_vec_len
% 			match_pt1(m,:)= [NaN, NaN];
% 		end
% 	end
% 	if(length(match_pt2)<max_vec_len)
% 		for m = length(match_pt2)+1 : max_vec_len
% 			match_pt2(m,:)= [NaN, NaN];
% 		end
% 	end
% 	match_pt1_tensor = cat(3,match_pt1_tensor,match_pt1);
% 	match_pt2_tensor = cat(3,match_pt2_tensor,match_pt2);
% end

% for image_ = 1:num_images-1
% 	match_pt1 = [];
% 	match_pt2 = [];
% 	for i = 1:length(match_pt1_tensor(:,1,1));
% 		match_pt1(i,:) = match_pt1_tensor(i,:,image_);
% 		match_pt2(i,:) = match_pt2_tensor(i,:,image_);
% 	end
% 	match_pt1 = match_pt1(~any(isnan(match_pt1), 2),:);
% 	% swap = match_pt1(:,1);
% 	% match_pt1(:,1) = match_pt1(:,2);
% 	% match_pt1(:,2) = swap;
% 	match_pt2 = match_pt2(~any(isnan(match_pt2), 2),:);
% 	% swap = match_pt2(:,1);
% 	% match_pt2(:,1) = match_pt2(:,2);
% 	% match_pt2(:,2) = swap;

% 	figure
% 	idk = showMatchedFeatures(image_tensor(:,:,image_),image_tensor(:,:,image_+1),match_pt1, match_pt2,'montage');
% % 	title('matched features');
% % end


% %gauss = fspecial('gaussian', window_size, sigma);

% function dist = cal_Euclidean_Distance(x,y,a,b)
% 	dist = (x-a)^2+(y-b)^2;
end

