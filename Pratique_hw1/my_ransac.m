%%% RANSAC

%Ransac parameters
ransac_iterations = 200;		%number of iteration
ransac_threshold = 3; 		%threshold
ransac_ratio = 0.6;			%ratio of inliers required to assert
							%that a model fits well to data
					


%loading given data set into variables
data_1 = load('hw1_data/data1.mat');
data_2 = load('hw1_data/data2.mat');
data_3 = load('hw1_data/data3.mat');

n_samples = max(size(data_1.pts));		%number of input points

%splitting data_1 into vectors
X_1 = data_1.pts(1,:).';
Y_1 = data_1.pts(2,:).';
%splitting data_2 into vectors
X_2 = data_2.pts(1,:).';
Y_2 = data_2.pts(2,:).';
%splitting data_3 into vectors
X_3 = data_3.pts(1,:).';
Y_3 = data_3.pts(2,:).';


%iterative model
ratio = 0;
model_m = 0;
model_c = 0;

%perform RANSAC iteration
for i = 1:ransac_iterations
	%pick any two random points
	n = 2;

	indice_1 = randi(n_samples,1)
	indice_2 = randi(n_samples,1)
	k = data_1.pts(:,indice_1)
	l = data_1.pts(:,indice_2)
	maybe_points = [data_1.pts(:,indice_1) data_1.pts(:,indice_2)];
	[slope, c] = model_of_line(maybe_points);
	x_y_inliers = [];
	num = 0;

	%find orthogonal lines to the model for all testing points
	for j = 1:n_samples
		if (j ~=  indice_1 ||  j ~=  indice_1)
			x0 = data_1.pts(1,j);
			y0 = data_1.pts(2,j);

			%find an intercept point of the model with a normal from point (x0,y0)
			[x1 y1] = find_intercept_point(slope, c, x0, y0);

			% distance from point to the model
        	perpendicular_dist = sqrt((x1 - x0)^2 + (y1 - y0)^2);

        	%check whether it's an inlier or not
        	if perpendicular_dist < ransac_threshold
        	    ponits_to_vector = [x0;y0];
            	x_y_inliers = [x_y_inliers ponits_to_vector];
            	num = num + 1;
        	end
        end
    end

    % in case a new model is better - save parameters
    if num/n_samples > ratio
    	ratio = num/n_samples;
    	model_m = slope;
    	model_c = c;
	end
	%X = sprintf('%s will be %d this year.',name,age);
	%disp(X)
	sprintf('Inlier ratio = %d',num/n_samples);
	sprintf('model_m = %d',model_m);
	sprintf('model_c = %d',model_c);

	%plot the current step
	plot_ransac(model_m,model_c,data_1);
	if num > n_samples*ransac_ratio
		sprintf('The model is found');		
		break;
	end
end

plot_ransac(model_m,model_c,data_1);
% scatter(data_1.pts(1,:),data_1.pts(2,:));
% hold on;
% x_plot_min = min(data_1.pts(1,:));
% x_plot_max = max(data_1.pts(1,:));
% x_plot = x_plot_min:x_plot_max;
% y_plot = model_m*x_plot+model_c;
% plot(x,y);
% hold off;


function [x y] = find_intercept_point(m, c, x0, y0)
	%{
		find the perpendicular distance of the point(x0,y0)
		from the selceted line model.
		find an intercept point of the line model with
    	a normal from point (x0,y0) to it
	    parameter m slope of the line model
	    parameter c y-intercept of the line model
	    parameter x0 point's x coordinate
	    parameter y0 point's y coordinate
	    output intercept point
	%}
	x = (x0 + m*y0 -m*c)/(1 + m^2);
	y = (m*x0 + (m^2)*y0 - (m^2)*c)/(1 + m^2) + c;
end


%Selecting random minimum points from datasets
function [slope,c] = model_of_line(points)
	% find a line model for the given points
	% parameter 'points' are slected points in vectorized form for model fitting
	% output = slope of line and constant c
	% ***Caveat***
	% vertical and horizontal lines should be treated differently
    % here we just add some noise to avoid division by zero
    
    % slope (gradient) of line
    slope = (points(2,2)-points(2,1)) / (points(1,2)-points(1,1) + eps); %added a little bit of noise to avoid division by zero or infinity
	c = points(2,2) - slope*points(2,1);
	
end

function plot_ransac(m,c,data)
	scatter(data.pts(1,:),data.pts(2,:));
	hold on;
	x_plot_min = min(data.pts(1,:));
	x_plot_max = max(data.pts(1,:));
	x_plot = x_plot_min:x_plot_max;
	y_plot = m*x_plot+c;
	plot(x_plot,y_plot);
	hold off;
end