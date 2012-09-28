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
    for j=2:1:length(x)-1
        px1 = diff(x(j-1),2);
        px2 = diff(x(j),2);
        px3 = diff(x(j+1),2);
        
        py1 = diff(y(j-1),2);
        py2 = diff(y(j),2);
        py3 = diff(y(j+1),2);
        
        %elastic(j-1).x = a*(px1-2.*px2+px3);
        %elastic(j-1).y = a*(py1-2.*py2+py3);
        
        elastic(j-1).x = 0;
        elastic(j-1).y = 0;
    end
    
    v = eval(U_times_B);
    
    Bm = diag([v(2),v(2),v(2)]);
    Bm(1,2) = v(3);
    Bm(2,3) = v(3);
    Bm(2,1) = v(1);
    Bm(3,2) = v(1);
    
    invBm= Bm^(-1);
    
    T = struct([]);
    t = 1;
    force = struct([]);
    for i=1:1:length(x)-2
        force(i).x = zi(i).zi + elastic(i).x;
        force(i).y = zi(i).zi + elastic(i).y;
        T(i).x = t.*invBm.*force(i).x;
        T(i).y = t.*invBm.*force(i).y;
    end   
end
m = img;