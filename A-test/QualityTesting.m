clc
close all
clear all

orginalStarry = imread('Images/starryNight.jpg');
starryBiggest = imread('Tests/starryNight20colorsLeft.png');
strarrySmallest = imread('Tests/starryNight20colorsLeftSmall.png');

orginalSmall = imresize(orginalStarry, 0.7);

[ssimval,ssimmap] = ssim(starryBiggest,orginalSmall);
imshow(ssimmap,[])
title(['Local SSIM Map with Global SSIM Value: ',num2str(ssimval)])

[ssimval,ssimmap] = ssim(strarrySmallest,orginalSmall);
figure;
imshow(ssimmap,[])
title(['Local SSIM Map with Global SSIM Value: ',num2str(ssimval)])

%% Kod för att lägga till kontur på orginalbilden
% a = i;
% b = repmat(mask, 1, 1, 3);
% a(mod(b,2)==1) = 0;
% imshow(a);
