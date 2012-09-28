function [seg_img] = clusters(img, num_class)

%% Getting Image
img_dir = '../img';
addpath(img_dir);
filename = fullfile(img_dir,'football.jpeg');
img = imread(filename);
img = rgb2gray(img);
figure, imshow(img,[]);

seg_img = img;