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
for iteration=1:1:1
    B = 0.5*[1,1,0;-2 2 0;1 -2 1];
    
    u = sym('u');
    U = [1, u, u^2];
    
    U_times_B =  U * B;
    
    u2eval = 20;
    u = 0:1/u2eval:1-1/u2eval;
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
    [mx,my]  = meshgrid(1:483,1:336);
    
    u = 11;
    zi = struct([]);
    for i=1:1:length(x)-2
        zi(i).zi = interp2(mx,my,Fx,C_s(i).x(u),C_s(i).y(u));
        zi(i).x = C_s(i).x(u);
        zi(i).y = C_s(i).y(u);
    end
    
    a = 1;
    elastic = struct([]);
    ddx = gradient(gradient(x));
    ddy = gradient(gradient(y));
        
    for j=2:1:length(ddx)-1
        elastic(j-1).x = a*(ddx(j-1)-2.*ddx(j)+ddx(j+1));
        elastic(j-1).y = a*(ddy(j-1)-2.*ddy(j)+ddy(j+1));
    end
    
    v = eval(U_times_B);
    
    Bm = diag([v(2),v(2),v(2)]);
    Bm(1,2) = v(3);
    Bm(2,3) = v(3);
    Bm(2,1) = v(1);
    Bm(3,2) = v(1);
    
    invBm= Bm^(-1);
    
    t = 1;
    sum_x = zeros(1,length(x)-2);
    sum_y = zeros(1,length(x)-2);
    for i=1:1:length(x)-2
        sum_x(i) = zi(i).zi + elastic(i).x;
        sum_y(i) = zi(i).zi + elastic(i).y;
    end
    
    new_x = t*invBm*sum_x';
    new_y = t*invBm*sum_y';
    
    %x = new_x;
    %y = new_y;
end
m = img;