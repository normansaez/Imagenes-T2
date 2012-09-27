%function [m] = snake(img)
clear all
%%
img_dir = '../img';
addpath(img_dir);
filename = fullfile(img_dir,'football.jpeg');
img = imread(filename);
img = rgb2gray(img);
%%
imshow(img,[]);
[x,y] = getpts;
% XXX: Al menos 3 puntos de control por B-spline

B = 0.5*[1,1,0;-2 2 0;1 -2 1];
u = sym('u');
U = [u^2, u, 1];

U_times_B =  U * B;

u2eval = 20;
u = 0:1/u2eval:1-1/u2eval;

C_u_vals = struct([]);
for j=1:1:length(x)-2
    
    px(1) = x(j);
    px(2) = x(j+1);
    px(3) = x(j+2);
    
    py(1) = y(j);
    py(2) = y(j+1);
    py(3) = y(j+2);
    
    Cx =  U_times_B * px';
    Cy =  U_times_B * py';    
    
    C_u_vals(j).x = eval(Cx);
    C_u_vals(j).y = eval(Cy);
    
end
m = img;
