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
end
