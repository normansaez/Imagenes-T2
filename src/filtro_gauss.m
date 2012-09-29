function [img] = filtro_gauss(imagen)
h = fspecial('gaussian', [3 3], 0.5);
img = imfilter(imagen,h);