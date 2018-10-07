path = dir('hw2_imgs/*.jpg');
nfiles = length(path);

window_size = 3;
corner_threshold = 7000000000;

for i = 1:nfiles
    currImagePath=fullfile(path(i).folder, path(i).name);
    %read the image
    image_=imread(currImagePath);
    myharris(image_,window_size,corner_threshold);
end