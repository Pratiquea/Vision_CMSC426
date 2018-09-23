
path=dir('train_images/*.jpg');
nFiles=length(path);

%vert stack to store R,G,B channels of entire dataset

orange_pixels_lab = [];	%zeros(1,3);
% imageStack_lab = [];%zeros(1,3);
for i=1:nFiles-1
    %currImagePath=path(i).folder+path(i).name;
    currImagePath=fullfile(path(i).folder, path(i).name);
    %read the image
    currImage=imread(currImagePath);
    %converting image from RGB to l*a*b*   
    currImage_lab = rgb2lab(currImage);
    % Define thresholds for channel l for lab image
    channel_l_Min = 19.510;
    channel_l_Max = 91.394;

    % Define thresholds for channel a for lab image
    channel_a_Min = 14.923;
    channel_a_Max = 51.449;
    % Define thresholds for channel b for lab image
    channel_b_Min = 10.049;
    channel_b_Max = 46.717;

    % Create mask based on chosen histogram thresholds
    %(currImage_lab(:,:,1) >= channel_l_Min ) & (currImage_lab(:,:,1) <= channel_l_Max) & ...
    sliderBW_lab = (currImage_lab(:,:,2) >= channel_a_Min ) & (currImage_lab(:,:,2) <= channel_a_Max) & ...
                    (currImage_lab(:,:,3) >= channel_b_Min ) & (currImage_lab(:,:,3) <= channel_b_Max);
    BW_lab = sliderBW_lab;

    % Initialize output masked image based on input image.
    maskedRGBImage_lab = currImage_lab;
    % Set background pixels where BW is false to zero.
    maskedRGBImage_lab(repmat(~BW_lab,[1 1 3])) = 0;
    imshow(maskedRGBImage_lab);
    maskedRGBImage_lab=reshape(maskedRGBImage_lab,640*480,3);
    %gathering all orange pixels(lab) in an matrix
    for pixels = 1:(640*480)
        if BW_lab(pixels) == 1
            orange_pixels_lab = vertcat(orange_pixels_lab,maskedRGBImage_lab(pixels,:));
        end
    end
end


[mean_1,sigma_1] = gmm_train(orange_pixels_lab,5);


function [mu sigma] = gmm_train(x, k)
	% x is nx3 data set where each coloumn are RGB or lab or hsv
	% k is the number of clusters
	num_pts = size(x,1);
	dim = size(x,2);
	%initializing parameters
	error_convergence = 0.001;
	pi_ = randn(k,1);
	% initializing mean(k x dim matrix)
	mean_ = 200*(rand(k,dim));
	prev_mean = mean_;
	%initializing covariance(dim x dim x k matrix)
	%reference matlab answers
	covariance = 100*(reshape(repmat(diag(ones(dim,1)),1,k),[dim, dim, k]));
	%threshold for comparing means
	threshold = 1;
	if max(max(x))<10
    threshold = 10e-5;
	end
	%posteriors(n x k matrix)
	posteriors = zeros(num_pts,k);
	%maximum training iterations
	max_iters = 10;

	iters = 0;
	while iters<max_iters
		%E-step
		for i = 1:k
			%posteriors = z x n matrix
			posteriors(:,i) = pi_(i)*single_gaussian_predict_pixel(mean_(i,:),covariance(:,:,k),x);
		end
		%summing the posteriors of all clusters : nx1
		sum_posteriors = sum(sum(posteriors(:,i),2),1);
		%normalizing the posterior vector(again n x 1 matrix)
		posteriors = posteriors/sum_posteriors;

		%M-step
		for i = 1:k
			mean_(i,:) = sum((posteriors(:,i).*x))/sum(posteriors(:,i));
			cov_posterior_prod = zeros(3,3);
			posteriors_sum = 0;
			for j = 1:size(x,1)
				cov_posterior_prod =cov_posterior_prod + posteriors(j,i)*(((x(j,:)-mean_(i,:))')*(x(j,:)-mean_(i,:)));
				posteriors_sum = posteriors_sum + posteriors(j,i);
			end
			covariance(:,:,i) = cov_posterior_prod/posteriors_sum;
			pi_(i) = sum(posteriors(:,i))/size(x,1);
        end

		iters = iters+1;
	end

	norm(prev_mean - mean_)
	if norm(prev_mean - mean_) < threshold
		return
	end
	prev_mean = mean_;

    plot(x(:,1),x(:,2),x(:,3));
	hold on
	error_ellipse(covariance(:,:,1),mean_(1,:))
	for i = 2:k
		error_ellipse(covariance(:,:,i),mean_(i,:));
	end
	hold off

end 		%function end


%funtion to generate random PSD covariance matrix
%		**reference: Matlab answers
% function M=randCov(N)
% d = 10*rand(N,1); % The diagonal values
% t = triu(bsxfun(@min,d,d.').*rand(N),1); % The upper trianglar random values
% M = diag(d)+t+t.'; % Put them together in a symmetric matrix
% end


%funtion to generate random PSD covariance matrix
%		**reference: math stackexchange
% function A = generateSPDmatrix(n)
% % Generate a dense n x n symmetric, positive definite matrix
% A = rand(n,n); % generate a random n x n matrix
% % construct a symmetric matrix using either
% A = 0.5*(A+A'); OR
% % A = A*A';
% % The first is significantly faster: O(n^2) compared to O(n^3)
% % since A(i,j) < 1 by construction and a symmetric diagonally dominant matrix
% %   is symmetric positive definite, which can be ensured by adding nI
% A = A + n*eye(n);
% end


function [posteriors] = single_gaussian_predict_pixel(mean_, cov_, x)
	% x: rgb 1x3
	% cov = 
	% mu: 1x3
	dim = size(x, 2);
	% x_mu = bsxfun(@minus, x, mu);
	cov_orange_lab = cov_;
	mean_orange_lab = mean_;

	likelihood_orange_lab = mvnpdf([x(:,1) x(:,2) x(:,3)], mean_orange_lab, cov_orange_lab);
    posteriors = likelihood_orange_lab;
    % filtered_img = reshape(posterior_orange_lab,640,480);
ende