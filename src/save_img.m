function [] = save_img(img, dir, name)
imagename = fullfile(dir,sprintf('%s.eps',name));
f = figure; 
imshow(img,[]);
set(f,'name',imagename,'numbertitle','off')
print(gcf,'-dpsc2',imagename);
%close;