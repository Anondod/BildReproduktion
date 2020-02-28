% org = original image (in RGB)
% rep = reproduced image (in RGB)
% Obs! images needs to be of same dimensions

% resW, resH = screens resolution
% dist = viewingdistance
function [qual,avg] = sCIELabMetric(org,rep,resW,resH,screenSize,dist)

subplot(1,3,1)
imshow(org)
title('Original')
subplot(1,3,2)
imshow(rep)
title('Reproduced')

% Convert images to XYZ for scielab()
orgXYZ = rgb2xyz(org);
repXYZ = rgb2xyz(rep);

% Calculate screens ppi 
% and then samples per degree from viewing distance dist
ppi = sqrt(resW^2+resH^2)/screenSize;
sampPerDeg = ppi*dist*tan(pi/180);

% Get S-CIELab quality measure
% Vector = whitepoint for CIED65 in XYZ
qual = scielab(sampPerDeg, orgXYZ, repXYZ, [95.05, 100, 108.9], 'xyz');
avg = mean(mean(qual));

% Prints quality result. Obs! not normalized
subplot(1,3,3)
imshow(qual)
title('S-CIELab quality diff')

end

