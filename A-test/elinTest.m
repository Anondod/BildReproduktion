% expansion test
clc
close all

% 0 draws nothing
% 1 draws every 1% of total segments finished
% 2 draws when each segment is finished
% 3 shows every pixel update
drawMode = 2;

% have 1px contours around segments?
doContours = 1;

% Open/close structuring element (STREL) 
closeSphere = strel('sphere', 3);
closeLineHorizontal = strel('line', 50, 0);
closeLineVertical = strel('line', 50, 90);


% for step 2 and the total coverage at end
drawcol(1,1,:) = [0,0,1];

distCityBlock = @(p1,p2) sum(abs(p2-p1));
distEuclidean = @(p1,p2) sqrt(sum((p2-p1).^2));

% choose one from above: euclidean , cityblock
distFun = distEuclidean;
distanceFac = 0.4;

qSize = 4000;

maxErr = 8000;

% if a pixel of the segment has to large an error the segment stops
cancelThreshold = 80;

%nr mosaics (max)
n = 5000;

% DECIDE INPUT IMAGE HERE
I = im2double(imresize(imread('elin.jpg'),.8));

% Blur image?
I = imgaussfilt(I,2);
Isize = [size(I,1) size(I,2)];

Ilab = rgb2lab(I);

Res = zeros(size(I));

mask = uint16(zeros(Isize(1),Isize(2)));

q = zeros(3, qSize);


printmod = ceil(n/50);
printmod = printmod + mod(printmod,2);

covered = 0;


% main segment for-loop
for i=2:2:2*n

% alt faster ver. but looks worse maybe. also breaks coverage
%[py_temp,px_temp] = find(mask==0,1,'first');
[py_temp,px_temp] = find(mask==0);

unfilled = size(px_temp,1);

if(unfilled==0)
    break;
end

%progress bar
if(mod(i,printmod)==0)
    covered = 1-unfilled/(Isize(1)*Isize(2));
    fprintf(i/2+"/"+n+" ("+100*i/(2*n)+"%%) segments done\n");
    fprintf(100*(covered)+"%% coverage\n\n");
    
    if(drawMode == 1)
        warning off;
        imshow(Res)
        warning on;
    end
    drawnow
end

index = randi(unfilled);
py = py_temp(index);
px = px_temp(index);

color = I(py,px,:);
e = Ilab - Ilab(py,px,:);
E = sqrt(e(:,:,1).^2+e(:,:,2).^2+e(:,:,3).^2);

qUsed = 1;
q(:,1) = [py px E(py,px)];

drawValue = 0;
% Loop per segment
acc_error = 0;
while acc_error < maxErr
    [err,lowest] = min(q(3,1:qUsed));
    p = q(1:2,lowest);
    acc_error = acc_error + err;
    
    if(qUsed ~= 0 && err < cancelThreshold)
        q(:,lowest) = q(:,qUsed);
        qUsed = qUsed-1;
    else
        break;
    end
        
    mask(p(1),p(2)) = i;
    
    % adding adjacents to queue and mask
    
    % pdist substitute: (sqrt((px-pt(2)).^2+(py-pt(1)).^2))

    
    pt =[p(1)+1, p(2)];
    if(p(1) ~= Isize(1) && mask(pt(1),pt(2))==0)
        mask(pt(1),pt(2)) = i+1;
        qUsed = qUsed+1;
        q(:,qUsed) = [pt(1),pt(2), E(pt(1),pt(2))+distFun([pt(1),pt(2)],[py,px])*distanceFac];
    end
    pt =[p(1)-1, p(2)];
    if(p(1) ~= 1 && mask(pt(1),pt(2))==0)
        mask(pt(1),pt(2)) = i+1;
        qUsed = qUsed+1;
        q(:,qUsed) = [pt(1),pt(2), E(pt(1),pt(2))+distFun([pt(1),pt(2)],[py,px])*distanceFac];
    end
    pt =[p(1), p(2)+1];
    if(p(2) ~= Isize(2) && mask(pt(1),pt(2))==0)
        mask(pt(1),pt(2)) = i+1;
        qUsed = qUsed+1;
        q(:,qUsed) = [pt(1),pt(2), E(pt(1),pt(2))+distFun([pt(1),pt(2)],[py,px])*distanceFac];
    end
    pt =[p(1), p(2)-1];
    if(p(2) ~= 1 && mask(pt(1),pt(2))==0)
        mask(pt(1),pt(2)) = i+1;
        qUsed = qUsed+1;
        q(:,qUsed) = [pt(1),pt(2), E(pt(1),pt(2))+distFun([pt(1),pt(2)],[py,px])*distanceFac];
    end
   
    
    % doesnt look good when it draws around lower or right edge atm
    if(drawMode == 3 && mod(drawValue, 10)==0)
        width = 100;
        scope1 = max(1,p(1)-width);
        scope1 = scope1:min(Isize(1),scope1+2*width);
        scope2 = max(1,p(2)-width);
        scope2 = scope2:min(Isize(2),scope2+2*width);
        
        draw = I(scope1,scope2,:) + double(mask(scope1,scope2) == i).*drawcol;
        %draw = E(scope1,scope2,:)/100 + double(mask(scope1,scope2) == i).*drawcol;
        imshow(draw);
        
        %q(:,1:qUsed+1) % what's in q at the moment
        pause(0.001);
    end
    
    drawValue = drawValue+1;
end


if(doContours==0)
    mask(mask==i+1) = 0;
end

tempClose = mask == i;

tempClose = imclose(tempClose, closeLineHorizontal);
tempClose = imclose(tempClose, closeLineVertical);
tempClose = imclose(tempClose, closeSphere);

sphere = strel('sphere', 2);

tempDilate = imdilate(tempClose, sphere);

contorsResult = logical(tempDilate - tempClose);

mask(tempClose) = i;
mask(contorsResult) = i+1;




mask3d = repmat(mask,1,1,3);
Res(mask3d==i) = color(1,1,1);
mask3d(:,:,1)=0;
Res(mask3d==i) = color(1,1,2);
mask3d(:,:,2)=0;
Res(mask3d==i) = color(1,1,3);


Res(mask3d==i+1) = 1;
mask3d(:,:,1)=0;
Res(mask3d==i+1) = 1;
mask3d(:,:,2)=0;
Res(mask3d==i+1) = 1;

if(drawMode == 2)
    imshow(Res)
end

end

fprintf("done at " + 100*(1-unfilled/(Isize(1)*Isize(2)))+"%% coverage\n");

% imclose for better looking segments maybe
%a = imclose(mask==1,strel('disk',10));

%%
close all
figure;
draw = I + double(mask).*drawcol;
imshow(draw);
title("total coverage");

figure;
imshow(double(mask)/i)
title("segments");


figure;
imshow(Res)
title("result");


% SAVE COMMAND
%imwrite(Res,'insertname.png')