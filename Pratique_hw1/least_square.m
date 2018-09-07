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
%defining a vector of ones for convenience
% vector_of_ones = ones(max(size(X_1)),1);
% A = [X_1 Y_1 vector_of_ones];
%optimization parameter
% theta = [a;b;c];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%			Homework question 2 - Linear Regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%using (backslash)
format long
%for data_1
B = inv(X_1.'*X_1)*X_1.'*Y_1;  % can also be written as B = X_1\Y_1;
y_cal = B*X_1;
figure('Name','linear regression using backslash');
subplot(2,2,1);
scatter(X_1,Y_1);
hold on
plot(X_1,y_cal);
hold off
%for data_2
B = X_2\Y_2;
y_cal = B*X_2;
subplot(2,2,2);
scatter(X_2,Y_2);
hold on
plot(X_2,y_cal);
hold off
%for data_3
B = X_3\Y_3;
y_cal = B*X_3;
subplot(2,2,3);
scatter(X_3,Y_3);
hold on
plot(X_3,y_cal);
hold off




%using polyfit
p = polyfit(X_1,Y_1,1); 
f = polyval(p,X_1); 
figure('Name','using polyfit function');
subplot(2,2,1);
plot(X_1,Y_1,'o',X_1,f,'-');
p = polyfit(X_2,Y_2,1); 
f = polyval(p,X_2); 
subplot(2,2,2);
plot(X_2,Y_2,'o',X_2,f,'-');
p = polyfit(X_3,Y_3,1); 
f = polyval(p,X_3); 
subplot(2,2,3);
plot(X_3,Y_3,'o',X_3,f,'-');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%			Homework question 2 - Regularization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for data_1

% lambda = [-100000 -10000 -1000 -100 -10 1 10 100 1000 10000 100000];  %created a vector to reduce manual hit and trial effort
lambda = -100000;
% for i = 1:max(size(lambda))
	B = inv(X_1.'*X_1 + lambda)*X_1.'*Y_1; 
	% y_cal = B*X_1;
	y_cal_adjust = B*X_1+10;		%added an offset of ten for fun(fits better!)
	figure('Name',"Using Regularization");
	% subplot(2,1,1);
	% scatter(X_1,Y_1);
	% hold on
	% plot(X_1,y_cal);
	% hold off
	subplot(2,2,1);
	scatter(X_1,Y_1);
	hold on
	plot(X_1,y_cal_adjust);
	hold off
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for data_2
lambda = -280000;
B = inv(X_2.'*X_2 + lambda)*X_2.'*Y_2; 
y_cal_adjust = B*X_2+10;		%
subplot(2,2,2);
scatter(X_2,Y_2);

hold on
plot(X_2,y_cal_adjust);
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for data_3
lambda = -500000;
B = inv(X_3.'*X_3 + lambda)*X_3.'*Y_3; 
y_cal_adjust = B*X_3+10;		%
subplot(2,2,3);
scatter(X_3,Y_3);
hold on
plot(X_3,y_cal_adjust);
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%			Homework question 2 - RANSAC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for data_1