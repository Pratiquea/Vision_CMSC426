image_ = rgb2gray(imread('football.jpg'));
% subplot(1,2,1);
% subplot(1,2,2);
corner_score= cornermetric(image_);
% subplot(1,2,1);
% imagesc(corner_score);
% title('cornermetric');
local_maxima = imregionalmax(corner_score);
% figure
% subplot(1,2,2);
% imagesc(local_maxima);
% figure
% % subplot(1,2,1);
% surf(corner_score);
% figure
% im_irm = corner_score;
% im_irm(local_maxima) = 0;
% surf(im_irm);
N_best = 100;
N_strong= sum(local_maxima(:)==1);
[y x] = find(local_maxima);

figure;
imshow(image_);
hold on
plot(x,y,'ys');
title('original image with local maxima corners');
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Adaptive Non-Maxima Supression step %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% radial_dist = [];
radial_dist_sorted = [];
corner_score_sorting_param = [];
anms_output = [];
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
% [corner_score_sorting_param, sort_index] = sort(corner_score_sorting_param,'descend');
% radial_dist_sorted = sort(radial_dist,'descend');
cells = struct2cell(anms_output);
sortvals = cells(:,1,1);
mat = (cell2mat(sortvals));
[sorted_,ix] = sort(mat);
anms_output = anms_output(ix);
% [y1 x1] = find(corner_score_sorting_param(1:N_best));
% [y1 x1] = find(radial_dist_sorted(1:N_best));
x1 = anms_output(1:N_best).x;
y1 = anms_output(1:N_best).y;

figure;
imshow(image_);
hold on
plot(x1,y1,'ys');
title('original image after anms');
hold off

function dist = cal_Euclidean_Distance(x,y,a,b)
	dist = (x-a)^2+(y-b)^2;
end

