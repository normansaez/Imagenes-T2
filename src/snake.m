%function [m] = snake(img)

%% Getting Image
close all;
clear all;
img_dir = '../img';
addpath(img_dir);
filename = fullfile(img_dir,'football.jpeg');
img = imread(filename);
img = rgb2gray(img);
figure, imshow(img,[]);

imshow(img,[]);
[x,y] = getpts;

% Control poligon
x(length(x)+1) = x(1);
y(length(y)+1) = y(1);
x(length(x)+1) = x(2);
y(length(y)+1) = y(2);
% Control poligon

plot(x,y,'--r')
set(gca,'YDir','reverse');
hold on
%pause

B = 0.5*[1,1,0;-2 2 0;1 -2 1];

u = sym('u');
U = [1, u, u^2];

U_times_B =  U * B;

%u2eval = 20;
%u = 0:1/u2eval:1-1/u2eval;
u = linspace(0,1,20);

C_s = struct([]);
for j=2:1:length(x)-1
    
    px(1) = x(j-1);
    px(2) = x(j);
    px(3) = x(j+1);
    
    py(1) = y(j-1);
    py(2) = y(j);
    py(3) = y(j+1);
    
    
    Sx =  U_times_B * px';
    Sy =  U_times_B * py';

    %Number of segments depends on number of control points
    i = j - 1;
    fprintf('segment: %d\n',i)
    C_s(i).x = eval(Sx);
    C_s(i).y = eval(Sy);
    plot(C_s(i).x,C_s(i).y,'-b')
    set(gca,'YDir','reverse');
    hold on
    %pause
end

[Fx, Fy] = gradient(gradient(double(img)).^2);
m = img;
