mat_filename = 'mean_and_cov.mat';
%loading mean and covariance from mat file
data = load(mat_filename);

cov_orange_lab = data.cov_orange_lab;
mean_orange_lab = data.mean_orange_lab;
