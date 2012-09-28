function [seg_img] = clusters(img, k)

%% Getting Image
img_dir = '../img';
addpath(img_dir);
filename = fullfile(img_dir,'football.jpeg');
img = imread(filename);
img = rgb2gray(img);
figure, imshow(img,[]);

index = find(img(:));
centroids =(1:k)*max(img(:))/(k+1);
