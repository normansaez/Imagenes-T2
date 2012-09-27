%% Procesamiento Avanzado de Imagenes
%  Tarea N2
%  mailto: nfsaez@uc.cl
%%

%% Take init time and clear/close all vars
tic;
clear all;
close all;

%% Getting Image
img_dir = '../img';
addpath(img_dir);
filename = fullfile(img_dir,'football.jpeg');
img = imread(filename);
img = rgb2gray(img);
%figure, imshow(img,[]);

%% Pregunta 1
num_class = 1;
c_img = clusters(img,num_class);
save_img(img,img_dir,'clusters');
%% Pregunta 2
s_img = snake(img);
save_img(s_img,img_dir,'snake');

%% Pregunta 3
%balloon
b_img = balloon(img);
save_img(b_img,img_dir,'balloon');

%filtro_gauss
f_img = filtro_gauss(img);
save_img(f_img,img_dir,'filtro_gauss');

%% Take final time and print it.
totaltime = toc;
fprintf('\nExecution time %.2f[min] or %.2f [sec]\n', totaltime/60, totaltime);