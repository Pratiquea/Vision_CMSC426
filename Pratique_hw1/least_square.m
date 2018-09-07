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
scatter(X_1,Y_1);
axis([-150 150 -150 150]);
hold on
plot(X_1,y_cal);
hold off
%for data_2
B = X_2\Y_2;
y_cal = B*X_2;
subplot(3,1,2);
scatter(X_2,Y_2);
axis([-150 150 -150 150]);
hold on
plot(X_2,y_cal);
hold off
%for data_3
B = X_3\Y_3;
y_cal = B*X_3;
subplot(3,1,3);
scatter(X_3,Y_3);
axis([-150 150 -150 150]);
hold on
plot(X_3,y_cal);
hold off

%%%%%%%%using polyfit
p = polyfit(X_1,Y_1,1); 
f = polyval(p,X_1); 
figure('Name','using polyfit function');
subplot(3,1,1);
axis([-150 150 -150 150]);
plot(X_1,Y_1,'o',X_1,f,'-');
p = polyfit(X_2,Y_2,1); 
f = polyval(p,X_2); 
subplot(3,1,2);
axis([-150 150 -150 150]);
plot(X_2,Y_2,'o',X_2,f,'-');
p = polyfit(X_3,Y_3,1); 
f = polyval(p,X_3); 
subplot(3,1,3);
axis([-150 150 -150 150]);
plot(X_3,Y_3,'o',X_3,f,'-');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%			Homework question 2 - Regularization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% for data_1
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
	scatter(X_1(:,2:end),Y_1);
	hold on;
	plot(X_1(:,2:end),y_cal_adjust);
	hold off;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for data_2
	lambda = -300000;
	X_2 = [ones(size(X_2,1),1) X_2];

	B = inv(X_2.'*X_2 + lambda*eye(2))*X_2.'*Y_2; 
	y_cal_adjust = X_2*B;		%
	subplot(3,1,2);
	scatter(X_2(:,2:end),Y_2);
	axis([-150 150 -150 150]);
	hold on
	plot(X_2(:,2:end),y_cal_adjust);
	hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for data_3
	lambda = -500000;
	X_3 = [ones(size(X_3,1),1) X_3];

	B = inv(X_3.'*X_3 + lambda*eye(2))*X_3.'*Y_3; 
	y_cal_adjust = X_3*B;	
	subplot(3,1,3);
	scatter(X_3(:,2:end),Y_3);
	axis([-150 150 -150 150]);
	hold on
	plot(X_3(:,2:end),y_cal_adjust);
	hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%			Homework question 2 - RANSAC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for data_1