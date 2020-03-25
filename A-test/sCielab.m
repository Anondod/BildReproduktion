clc
clear all
orgImg = im2double(imread('cypress.JPG'));
repImg = im2double(imread('cypress20colors.png'));
subplot(1,3,1)
imshow(orgImg)
subplot(1,3,2)
imshow(repImg)

% Konvertera bilderna till XYZ för användning i funktionen
orgXYZ = rgb2xyz(orgImg);
repXYZ = rgb2xyz(repImg);

ppi = sqrt(1920^2+1200^2)/15; %kanske 20 inches, dålig på att uppskatta
sampPerDeg = ppi*20*tan(pi/180);

qual = scielab(sampPerDeg, orgXYZ, repXYZ, [95.05, 100, 108.9], 'xyz');
nearAvg = mean(mean(qual));
subplot(1,3,3)
imshow(qual)