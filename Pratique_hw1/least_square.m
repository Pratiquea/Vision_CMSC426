% loading given data set into variables
data_1 = load('hw1_data/data1.mat');
data_2 = load('hw1_data/data2.mat');
data_3 = load('hw1_data/data3.mat');
%splitting data_1 into vectors
X_1 = data_1.pts(1,:).';
Y_1 = data_1.pts(2,:).';
%splitting data_2 into vectors
X_2 = data_2.pts(1,:).';
Y_2 = data_2.pts(2,:).';
%splitting data_3 into vectors
X_3 = data_3.pts(1,:).';
Y_3 = data_3.pts(2,:).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%			Homework question 2 - Linear Regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%using (backslash)
format long
%for data_1
B = inv(X_1.'*X_1)*X_1.'*Y_1;  % can also be written as B = X_1\Y_1;
y_cal = B*X_1;
figure('Name','linear regression using backslash');
subplot(3,1,1);
plot(X_1,Y_1,'o');
axis([-150 150 -150 150]);
hold on
plot(X_1,y_cal,'LineWidth',2);
hold off
%for data_2
B = X_2\Y_2;
y_cal = B*X_2;
subplot(3,1,2);
plot(X_2,Y_2,'o');
axis([-150 150 -150 150]);
hold on
plot(X_2,y_cal,'LineWidth',2);
hold off
%for data_3
B = X_3\Y_3;
y_cal = B*X_3;
subplot(3,1,3);
plot(X_3,Y_3,'o');
axis([-150 150 -150 150]);
hold on
plot(X_3,y_cal,'LineWidth',2);
hold off

%%%%%%%%using polyfit
p = polyfit(X_1,Y_1,1); 
f = polyval(p,X_1); 
figure('Name','using polyfit function');
subplot(3,1,1);
axis([-150 150 -150 150]);
plot(X_1,Y_1,'o');
hold on;
plot(X_1,f,'-','LineWidth',2);
hold off;
p = polyfit(X_2,Y_2,1); 
f = polyval(p,X_2); 
subplot(3,1,2);
axis([-150 150 -150 150]);
plot(X_2,Y_2,'o')
hold on;
plot(X_2,f,'-','LineWidth',2);
hold off;
p = polyfit(X_3,Y_3,1); 
f = polyval(p,X_3); 
subplot(3,1,3);
axis([-150 150 -150 150]);
plot(X_3,Y_3,'o');
hold on;
plot(X_3,f,'-','LineWidth',2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%			Homework question 2 - Regularization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% for data_1
% lambda = [-1000000 -100000 -10000 1000 10000 100000];  %created a vector to reduce manual hit and trial effort
lambda = -100000;
X_1 = [ones(size(X_1,1),1) X_1];
% for i = 1:max(size(lambda))
	% disp(size(X_1));
	% figure(i);
	B = inv(X_1.'*X_1 + lambda*eye(2))*X_1.'*Y_1;
	y_cal_adjust = X_1*B;						% y = mx + c
	figure('Name',"Using Regularization");
	subplot(3,1,1);
	axis([-150 150 -150 150]);
	plot(X_1(:,2:end),Y_1,'o');
	hold on;
	plot(X_1(:,2:end),y_cal_adjust,'LineWidth',2);
	hold off;
% end
%%%%%%%%%%%%%%%%%%%% for data_2
	lambda = -300000;
	X_2 = [ones(size(X_2,1),1) X_2];

	B = inv(X_2.'*X_2 + lambda*eye(2))*X_2.'*Y_2; 
	y_cal_adjust = X_2*B;		%
	subplot(3,1,2);
	plot(X_2(:,2:end),Y_2,'o');
	axis([-150 150 -150 150]);
	hold on
	plot(X_2(:,2:end),y_cal_adjust,'LineWidth',2);
	hold off
%%%%%%%%%%%%%%%%%%%% for data_3
	lambda = -500000;
	X_3 = [ones(size(X_3,1),1) X_3];

	B = inv(X_3.'*X_3 + lambda*eye(2))*X_3.'*Y_3; 
	y_cal_adjust = X_3*B;	
	subplot(3,1,3);
	plot(X_3(:,2:end),Y_3,'o');
	axis([-150 150 -150 150]);
	hold on
	plot(X_3(:,2:end),y_cal_adjust,'LineWidth',2);
	hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%			Homework question 2 - RANSAC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'Line fitting using RANSAC');
%%%%%%%%%% for data_1
[model_m model_c inlierz] = my_ransacc(data_1, 100, 3, 90, false);
subplot(3,1,1);
plot_ransac(model_m,model_c, data_1);
title('Plot for data_1');

%%%%%%%%%% for data_2
[model_m model_c inlierz] = my_ransacc(data_2, 150, 4, 90,false);
subplot(3,1,2);
plot_ransac(model_m,model_c, data_2);
title('Plot for data_2');

%%%%%%%%%% for data_3
[model_m model_c inlierz] = my_ransacc(data_3, 200, 6, 90,false);
subplot(3,1,3);
plot_ransac(model_m,model_c, data_3);
title('Plot for data_3');

% Main ransac function 
function [model_slope, model_intercept, inliers] = my_ransacc(data, max_iterations, sigma_threshold, percentage_inliers, Visualize_RANCAC_plot)
	%Ransac parameters
	switch nargin
		case 1
			max_iterations = 200;			%number of iteration
			sigma_threshold = 3; 			%threshold
			percentage_inliers = 60;		%percentage(converted to ratio) of inliers to fit the model
			Visualize_RANCAC_plot = false;	%param to visually see RANSAC plot model at every iteration		
		case 2
			sigma_threshold = 3; 		%threshold
			percentage_inliers = 60;			%percentage(converted to ratio) of inliers to fit the model
			Visualize_RANCAC_plot = false;	%param to visually see RANSAC plot model at every iteration
		case 3
			percentage_inliers = 60;			%percentage(converted to ratio) of inliers to fit the model
			Visualize_RANCAC_plot = false;	%param to visually see RANSAC plot model at every iteration
		case 4
			Visualize_RANCAC_plot = false;	%param to visually see RANSAC plot model at every iteration
			
		case 5
			sprintf('received all 5 params');

		case num2cell(6:10)
			fprintf(['Kuch Zyada nhi ho gya?']);
		
		otherwise
			fprintf(['No parameters are provied\n', ...
					 'function call is as follows:\n', ...
					 'my_ransac(data, max_iterations, sigma_threshold, percentage_inliers)']);
			return;
	end
	percentage_inliers = percentage_inliers/100;
	% max_iterations = 200;		%number of iteration
	% sigma_threshold = 3; 		%threshold
	% percentage_inliers = 0.8;			%ratio of inliers required to assert that a model fits well to data
	n_samples = max(size(data.pts));		%number of input points
	window_open = false;

	%iterative model
	ratio = 0;
	model_slope = 0;
	model_intercept = 0;

	%perform RANSAC iteration
	for i = 1:max_iterations
		%pick any two random points
		n = 2;

		indice_1 = randi(n_samples,1);
		indice_2 = randi(n_samples,1);
		% k = data.pts(:,indice_1)
		% l = data.pts(:,indice_2)
		maybe_points = [data.pts(:,indice_1) data.pts(:,indice_2)];
		% calling a function to model the line using the randomly selected points
		[slope, c] = model_of_line(maybe_points);
		x_y_inliers = [];
		number_of_inliers = 0;

		%find orthogonal lines to the model for all testing points
		for j = 1:n_samples
			if (j ~=  indice_1 ||  j ~=  indice_2)
				x0 = data.pts(1,j);
				y0 = data.pts(2,j);

				%find an intercept point of the model with a normal from point (x0,y0)
				[x1 y1] = find_intercept_point(slope, c, x0, y0);

				% distance from point to the model
	        	perpendicular_dist = sqrt((x1 - x0)^2 + (y1 - y0)^2);

	        	%check whether it's an inlier or not
	        	if perpendicular_dist < sigma_threshold
	        	    ponits_to_vector = [x0;y0];
	            	x_y_inliers = [x_y_inliers ponits_to_vector];
	            	number_of_inliers = number_of_inliers + 1;
	        	end
	        end
	    end

	    % in case a new model is better - save parameters
	    if number_of_inliers/n_samples > ratio
	    	ratio = number_of_inliers/n_samples;
	    	model_slope = slope;
	    	model_intercept = c;
	    	inliers = number_of_inliers;
		end
		%X = sprintf('%s will be %d this year.',name,age);
		%disp(X)
		sprintf('Inlier ratio = %d',number_of_inliers/n_samples);
		sprintf('model_slope = %d',model_slope);
		sprintf('model_intercept = %d',model_intercept);

		if (Visualize_RANCAC_plot == true)
			if window_open ==true
				close 'Visualize RANSAC happening';
				window_open = false;
			end
			%plot the current step
			figure('Name','Visualize RANSAC happening');
			plot_ransac(model_slope,model_intercept,data);
			pause(0.05);
			window_open = true;
			% close 'Visualize RANSAC happening';

		end
		if number_of_inliers > n_samples*percentage_inliers
			sprintf('The model is found');		
			break;
		end
	end
end


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
	% ***Caveat**
	% vertical and horizontal lines should be treated differently
    % here we just add some noise to avoid division by zero
    
    % slope (gradient) of line
    slope = (points(2,2)-points(2,1)) / (points(1,2)-points(1,1) + eps); %added a little bit of noise to avoid division by zero or infinity
	c = points(2,2) - slope*points(2,1);
	
end

function plot_ransac(m,c,data)
	% scatter(data.pts(1,:),data.pts(2,:));
	plot(data.pts(1,:),data.pts(2,:),'o');
	hold on;
	axis([-150,150,-150,150]);
	x_plot_min = min(data.pts(1,:))*1.2;
	x_plot_max = max(data.pts(1,:))*1.2;
	x_plot = x_plot_min:x_plot_max;
	y_plot = m*x_plot+c;
	plot(x_plot,y_plot,'LineWidth',2);
	hold off;
end
