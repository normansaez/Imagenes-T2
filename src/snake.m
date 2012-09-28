%function [m] = snake(img)
%% Getting Image
close all;
clear all;
img_dir = '../img';
addpath(img_dir);
filename = fullfile(img_dir,'football.jpeg');
img = imread(filename);
img = rgb2gray(img);
%figure, imshow(img,[]);

imshow(img,[]);
[x,y] = getpts;

plot(x,y,'*r')
set(gca,'YDir','reverse');
hold on
pause

B = 0.5*[1,1,0;-2 2 0;1 -2 1];

u = sym('u');
U = [u^2, u, 1];

U_times_B =  U * B;

u2eval = 20;
u = 0:1/u2eval:1-1/u2eval;

C_s = struct([]);
for j=2:1:length(x)
        px(1) = x(j-1);
        px(2) = x(j);
        
        py(1) = y(j-1);
        py(2) = y(j);
        
    if j == length(x)
        px(3) = x(1);
        py(3) = y(1);
    else
        px(3) = x(j+1);
        py(3) = y(j+1);
    end
    
    Sx =  U_times_B * px';
    Sy =  U_times_B * py';
    
    i = j - 1;
    fprintf('Seg %d\n',i)
    C_s(i).x = eval(Sx);
    C_s(i).y = eval(Sy);
end

[Fx, Fy] = gradient(gradient(double(img)).^2);
m = img;
