%loading given data set into variables
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

%using linear regression (backslash)

format long
%for data_1
B = X_1\Y_1;
y_cal = B*X_1;
figure('Name','using linear regression (backslash)');
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
%for data_1
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