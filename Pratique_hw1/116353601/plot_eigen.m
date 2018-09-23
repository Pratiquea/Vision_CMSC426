%loading given data set into variables
data_1 = load('hw1_data/data1.mat');
data_2 = load('hw1_data/data2.mat');
data_3 = load('hw1_data/data3.mat');

%segregating data_1 into vectors
v1 = data_1.pts(1,:);
v2 = data_1.pts(2,:);
covarience_matrix_1 = cov(v1,v2);
[Evec_1, Eval_1] = eig(covarience_matrix_1)
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
plot(v1,v2,'o');
hold on;
x1 = [Evec_1(1,1),-Evec_1(1,1)]*sqrt(Eval_1(1,1)*2);	%scaling the eigen vectors by eigen values
y1 = [Evec_1(2,1),-Evec_1(2,1)]*sqrt(Eval_1(1,1)*2);

x2 = [Evec_1(1,2),-Evec_1(1,2)]*sqrt(Eval_1(2,2)*2);
y2 = [Evec_1(2,2),-Evec_1(2,2)]*sqrt(Eval_1(2,2)*2);
axis([-200,200,-200,200]);
plot(x1,y1,'LineWidth',3);
plot(x2,y2,'LineWidth',3);
hold off
% subplot(1,2,2);
% plot();

%scatter plotting the data_2
figure(2);
plot(v3,v4,'o');
hold on;
axis([-200,200,-200,200]);
x1 = [Evec_2(1,1),-Evec_2(1,1)]*sqrt(Eval_2(1,1)*2);	%scaling the eigen vectors by eigen values
y1 = [Evec_2(2,1),-Evec_2(2,1)]*sqrt(Eval_2(1,1)*2);

x2 = [Evec_2(1,2),-Evec_2(1,2)]*sqrt(Eval_2(2,2)*2);
y2 = [Evec_2(2,2),-Evec_2(2,2)]*sqrt(Eval_2(2,2)*2);
plot(x1,y1,'LineWidth',3);
plot(x2,y2,'LineWidth',3);
hold off;
%scatter plotting the data_3
figure(3);
plot(v5,v6,'o');
hold on;
axis([-200,200,-200,200]);
x1 = [Evec_3(1,1),-Evec_3(1,1)]*sqrt(Eval_3(1,1)*2);	%scaling the eigen vectors by eigen values
y1 = [Evec_3(2,1),-Evec_3(2,1)]*sqrt(Eval_3(1,1)*2);

x2 = [Evec_3(1,2),-Evec_3(1,2)]*sqrt(Eval_3(2,2)*2);
y2 = [Evec_3(2,2),-Evec_3(2,2)]*sqrt(Eval_3(2,2)*2);
plot(x1,y1,'LineWidth',3);
plot(x2,y2,'LineWidth',3);
hold off;


% tutorial in comments
% [u, s, v] = svd(a);
