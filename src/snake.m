%function [m] = snake(img)

%% Getting Image
close all;
clear all;
img_dir = '../img';
addpath(img_dir);
filename = fullfile(img_dir,'football.jpeg');
img = imread(filename);
img = rgb2gray(img);
figure(1), imshow(img,[]);

hold on
[x,y] = getpts;

% Control poligon
x(length(x)+1) = x(1);
y(length(y)+1) = y(1);
x(length(x)+1) = x(2);
y(length(y)+1) = y(2);
% Control poligon


for iteration=1:200
    figure(1)
    plot(x,y,'--r')
    set(gca,'YDir','reverse');
    hold on
    %pause
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
    if iteration == 1
        [Fx, Fy] = gradient(gradient(double(img)).^2);
        Fx = Fx/max(Fx(:));
        Fy = Fy/max(Fx(:));
    end
    [mx,my]  = meshgrid(1:483,1:336);
    
    u = 11;
    grad = struct([]);
    for i=1:1:length(x)-2
        grad(i).x = interp2(mx,my,Fx,C_s(i).x(u),C_s(i).y(u));
        grad(i).y = interp2(mx,my,Fy,C_s(i).x(u),C_s(i).y(u));
    end
    
    a = 1e0;
    elastic = struct([]);
    
    for j=2:1:length(x)-1
        elastic(j-1).x = a*(x(j-1)-2.*x(j)+x(j+1));
        elastic(j-1).y = a*(y(j-1)-2.*y(j)+y(j+1));
    end
    
    v = eval(U_times_B);
    
    Bm = diag([v(2),v(2),v(2)]);
    Bm(1,2) = v(3);
    Bm(2,3) = v(3);
    Bm(2,1) = v(1);
    Bm(3,2) = v(1);
    if iteration == 1
        invBm= Bm^(-1);
    end
    t = 0.1;
    sum_x = zeros(1,length(x)-2);
    sum_y = zeros(1,length(x)-2);
    for i=1:1:length(x)-2
        sum_x(i) = grad(i).x + elastic(i).x;
        sum_y(i) = grad(i).y + elastic(i).y;
    end
    
    
    
    for j=2:1:length(sum_x)-1
        sumx(1) = sum_x(j-1);
        sumx(2) = sum_x(j);
        sumx(3) = sum_x(j+1);
        
        sumy(1) = sum_y(j-1);
        sumy(2) = sum_y(j);
        sumy(3) = sum_y(j+1);
        
        delta_x = t*invBm*sumx';
        delta_y = t*invBm*sumy';
        
        %figure(3)
        %plot(iteration,delta_y)
        %hold on
        
        x(j-1) = x(j-1)   + delta_x(1);
        x(j)   = x(j)     + delta_x(2);
        x(j+1) = x(j+1)   + delta_x(3);
        
        y(j-1)   = y(j-1)   + delta_y(1);
        y(j)     = y(j)     + delta_y(2);
        y(j+1)   = y(j+1) + delta_y(3);
    end
end
figure, imshow(img,[])
hold on
plot(x,y,'--r')
hold on