%% functionname: function description
% function [outputs] = my_harris(arg)


i = imread('football.jpg');
figure(1);
subplot(2,3,1);
imshow(i);

% figure(2);
i2 = imgaussfilt(i,2);
i3 = imgaussfilt(i,3);
i4 = imgaussfilt(i,4);
subplot(2,3,2);
imshow(i3);
title('gauss 3');

BW2 = edge(rgb2gray(i2),'Canny');
subplot(2,3,3);
imshow(BW2);
title('gauss 2 Canny');

BW3 = edge(rgb2gray(i3),'Canny');
subplot(2,3,4);
imshow(BW3);
title('gauss 3 Canny');

BW4 = edge(rgb2gray(i4),'Canny');
subplot(2,3,5);
imshow(BW4);
title('gauss 4 Canny');

figure(2);
a = [-1 0 1; -1 4 1; -1 0 1];
b = imfilter(i,a,'conv');
imshow(b);
% end
