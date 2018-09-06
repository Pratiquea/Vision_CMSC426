%loading given data set into variables
data_1 = load('hw1_data/data1.mat');
data_2 = load('hw1_data/data2.mat');
data_3 = load('hw1_data/data3.mat');

% corrplot(data_1.pts);
% coef = pca(data_1.pts);
% plot(coef); 
%segregating data_1 into vectors
v1 = data_1.pts(1,:);
v2 = data_1.pts(2,:);
covarience_matrix_1 = cov(v1,v2);
[Evec_1, Eval_1] = eig(covarience_matrix_1)
visualize_eigen_1 = Evec_1(:,1)*(sqrt(Eval_1(1:1)));
visualize_eigen_2 = Evec_1(:,2)*(sqrt(Eval_1(4:4)));
%segregating data_1 into vectors
v3 = data_2.pts(1,:);
v4 = data_2.pts(2,:);
covarience_matrix_2 = cov(v3,v4);
[Evec_2, Eval_2] = eig(covarience_matrix_2);
%segregating data_1 into vectors
v5 = data_3.pts(1,:);
v6 = data_3.pts(2,:);
covarience_matrix_3 = cov(v5,v6);
[Evec_3, Eval_3] = eig(covarience_matrix_3);


%scatter plotting the data_1
figure(1);
% subplot(1,2,1);
scatter(v1,v2);
hold on;
% Evec_1=Evec_1.*100;
plot(visualize_eigen_1);
plot(visualize_eigen_2);
% quiver(Evec_1(1,1),Evec_1(1,2),230);
% quiver(1,1,230);
hold off
% subplot(1,2,2);
% plot();

%scatter plotting the data_2
figure(2);
scatter(v3,v4);
hold on;
plot(Evec_2);
hold off;
%scatter plotting the data_3
figure(3);
scatter(v5,v6);
hold on;
plot(Evec_3);
hold off;


% tutorial in comments
% [u, s, v] = svd(a);
