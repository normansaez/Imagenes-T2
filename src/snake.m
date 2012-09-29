%% Getting Image
close all;
clear all;
img_dir = '../img';
addpath(img_dir);
filename = fullfile(img_dir,'football.jpeg');
img = imread(filename);
img = rgb2gray(img);
figure(1)
imshow(img,[]);
hold on
[x,y] = getpts;

init_x = x;
init_y = y;
 
    
% Control poligon
x(length(x)+1) = x(1);
y(length(y)+1) = y(1);
x(length(x)+1) = x(2);
y(length(y)+1) = y(2);
% Control poligon

B = 0.5*[1,1,0;-2 2 0;1 -2 1];

u = sym('u');
U = [1, u, u^2];

U_times_B =  U * B;

u2eval = 20;
u = 0:1/u2eval:1-1/u2eval;
C_s = struct([]);
[Fx, Fy] = gradient(abs(gradient(double(img)).^2));
[mx,my]  = meshgrid(1:483,1:336);
u = 11;
v = eval(U_times_B);
u = 0:1/u2eval:1-1/u2eval;
Bm = diag([v(2),v(2),v(2)]);
Bm(1,2) = v(3);
Bm(2,3) = v(3);
Bm(1,3) = v(1);
Bm(2,1) = v(1);
Bm(3,2) = v(1);
Bm(3,1) = v(3);

invBm= Bm^(-1);

elasticidad = 1e-3;
viscosidad = 1e-2;
    
number_of_points = length(x);
number_of_segments = number_of_points-2;

for iteration=1:200
    figure(1)
    plot(x,y,'--r')
    set(gca,'YDir','reverse');
    hold on
    %pause
    
    for j=2:1:number_of_points-1
        
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
        plot(C_s(i).x,C_s(i).y,'-y')
        set(gca,'YDir','reverse');
        hold on
        %pause
    end
    u_index = 11;
    grad = struct([]);
    for i=1:1:number_of_segments
        grad(i).x = interp2(mx,my,Fx,C_s(i).x(u_index),C_s(i).y(u_index));
        grad(i).y = interp2(mx,my,Fy,C_s(i).x(u_index),C_s(i).y(u_index));
    end
    

    elastic = struct([]);
    
    for j=2:1:number_of_points-1
        elastic(j-1).x = elasticidad*(x(j-1)-2.*x(j)+x(j+1));
        elastic(j-1).y = elasticidad*(y(j-1)-2.*y(j)+y(j+1));
    end
    
    
    
    sum_x = zeros(1,number_of_segments);
    sum_y = zeros(1,number_of_segments);
    
    for i=1:1:number_of_segments
        sum_x(i) = grad(i).x + elastic(i).x;
        sum_y(i) = grad(i).y + elastic(i).y;
    end
    
    
    
    for j=2:1:number_of_segments-1
        sumx(1) = sum_x(j-1);
        sumx(2) = sum_x(j);
        sumx(3) = sum_x(j+1);
        
        sumy(1) = sum_y(j-1);
        sumy(2) = sum_y(j);
        sumy(3) = sum_y(j+1);
        
        delta_x = viscosidad*invBm*sumx';
        delta_y = viscosidad*invBm*sumy';
        
        %figure(3)
        %plot(iteration,delta_y)
        %hold on
        
        x(j-1) = x(j-1)   + delta_x(1);
        x(j)   = x(j)     + delta_x(2);
        x(j+1) = x(j+1)   + delta_x(3);
        
        y(j-1)   = y(j-1)   + delta_y(1);
        y(j)     = y(j)     + delta_y(2);
        y(j+1)   = y(j+1)   + delta_y(3);
    end
end
figure, imshow(img,[])
hold on
plot(x,y,'--r')
set(gca,'YDir','reverse');
hold on
for i=1:1:number_of_segments
    plot(C_s(i).x,C_s(i).y,'-y')
    set(gca,'YDir','reverse');
    hold on
end
plot(init_x,init_y,'*b')
set(gca,'YDir','reverse');
hold on
