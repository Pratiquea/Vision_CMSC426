%class tut
function out_img = apply_gaussian(img, mu, sig)

gauss = zeros(5);
for i = 1:4
	for j = 1:4
		x=i;y=j;
		gauss(i,j)= (1/(2*sig*sqrt(2*pi)) * exp(-1/2*( ((x-mu)^2 + (y-mu)^2)/2*sig^2))
	end
end
