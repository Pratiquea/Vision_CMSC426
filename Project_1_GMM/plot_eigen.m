%{
flow:
1. Gather the data in the form of 2*n matrix (2D data)
2. Get the covariance matrix
3. Eigenvectors and eigenvalues of this matrix will then be the basis for
the ellipse that forms.
4. The eigenvalues represent the magnitude of the major and minor axes and 
eigenvector the orientation.
%}

    
%load data
data1=load('data/data1.mat');
data2=load('data/data2.mat');
data3=load('data/data3.mat');

plotCov(data1.pts);
plotCov(data2.pts);
plotCov(data3.pts);

function plotCov(pts)

rawData=pts;
rawDatax=pts(1,:);
meanx=mean(rawDatax);
rawDatay=pts(2,:);
meany=mean(rawDatay);

%variance within and between variables - covariance
covData=cov(rawDatax,rawDatay);

%eigendecomposition of this matrix to get rotation and scaling matrix
 %s = -2 * log(1 - 0.9);
[U,V,F]=svd(covData);

%creating a stream of axes 
x=linspace(-1,1)
y=linspace(-1,1)

mean1=[meanx;meany];
D=[[x;zeros(size(x))],[zeros(size(y));y]]
u=(U * sqrt(V))*[D(1,:);D(2,:)];
scatter(rawDatax,rawDatay);
hold on

%plot(majAxis(1,:),majAxis(2,:),'xr')
%plot(minAxis,'xb')
plot(D(1,:),D(2,:),'xg')
plot(meanx+u(1,:),meany+u(2,:),'xb');
hold off
pause(3)
close
end
