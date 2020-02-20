% from https://www.projectrhea.org/rhea/index.php/Implementation_of_Color_Quantization

% uses k-means -> high amount of colors = very computing heavy

%%
clear all
close all
clc

% Amount of clusteres/colors
clusteres = 5;

% original image to color quantize
org = im2double(imread('elin.jpg'));

subplot(1,2,1)
imshow(org)
title('original')

% convert from RGB to work in CIELab
lab = rgb2lab(org);

% dela upp i olika färgkanaler
L = lab(:,:,1);
a = lab(:,:,2);
b = lab(:,:,3);

% reshape color-channels into column vectors for kmeans()
img(:,1) = L(:);
img(:,2) = a(:);
img(:,3) = b(:);

% perform kmeans-clustering
[clusterInd, centroidLocations] = kmeans(img, clusteres);

% find pairwise distances between original image and the quantized one
imgQ = pdist2(img, centroidLocations);

% find shortest distances
% minQ, min distance for every row
[minQ, nearestClusterInd] = min(imgQ,[],2);
% find assigned clusters' colors for every pixel
colQuant = centroidLocations(nearestClusterInd,:);

% reshape from vectors into image matrix again 
% to same size as the original image
[rows,cols,channels] = size(org);
result(:,:,1) = reshape(colQuant(:,1),rows,cols);
result(:,:,2) = reshape(colQuant(:,2),rows,cols);
result(:,:,3) = reshape(colQuant(:,3),rows,cols);

% convert back into RGB
result = lab2rgb(result);

subplot(1,2,2)
imshow(result)
title('color quantization')