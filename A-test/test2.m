clear all;

I = im2double(imread('test6.jpg'));

[grad_m,grad_dir] = imgradient(mean(I,3));
blur_m = imgaussfilt(grad_m, 1);

imshow(blur_m)

[px,py] = find(blur_m == min(min(blur_m)));


val = blur_m(py(1),px(1));
e = 0.1;
upper = val + e;
lower = val - e;

mask = blur_m > lower & blur_m < upper;

marker = false(size(mask));
marker(py(1),px(1)) = true;
reconImage = imreconstruct(marker, mask, 4);

%see the origin point maybe
%mask(py-5:py+5,px-5:px+5) = 0;

figure;
imshow(reconImage)

se = strel('disk',2);
morphed = imclose(reconImage,se);
figure;
imshow(morphed)

%%
blur = 4;

[grad_m,grad_dir] = imgradient(mean(I,3));
blur_m = imgaussfilt(grad_m, blur);
[grad_m,grad_dir] = imgradient(blur_m);
blur_m = imgaussfilt(grad_m, blur);

imshow(imgradient(blur_m))

%%
panes = zeros(size(mask));
n = 100;
x = floor(gallery('uniformdata',[1 n],0).*size(mask,2));
y = floor(gallery('uniformdata',[1 n],5).*size(mask,1));
m = voronoi2mask(x,y, size(mask));

imshow(m/n)
figure;
a = m~=0;
imshow(a)
