% expansion test
clear all;

% 0 draws nothing
% 1 draws when each segment is finished
% 2 shows every pixel update
drawMode = 0;

% have 1px contours around segments?
doContours = 1;

% for step 2 and the total coverage at end
drawcol(1,1,:) = [0,0,1];

distanceFac = 0.2;
% euclidean , cityblock
distMethod = 'euclidean';

qSize = 4000;

maxErr = 10000;

% if a pixel of the segment has to large an error the segment stops
cancelThreshold = 60;

%nr mosaics (max)
n = 5000;

% DECIDE INPUT IMAGE HERE
I = im2double(imresize(imread('test1.jpg'),.4));

% Blur image?
% I = imgaussfilt(I,1);
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

% alt faster ver. but looks worse maybe
%[py,px] = find(mask==0,1000,'first');
[py,px] = find(mask==0);

unfilled = size(px,1);

if(unfilled==0)
    break;
end

%progress bar
if(mod(i,printmod)==0)
    covered = 1-unfilled/(Isize(1)*Isize(2));
    fprintf(i/2+"/"+n+" ("+100*i/(2*n)+"%%) segments done\n");
    fprintf(100*(covered)+"%% coverage\n\n");
    drawnow
end

index = randi(unfilled);
py = py(index);
px = px(index);

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
    
    pt =[p(1)+1, p(2)];
    if(p(1) ~= Isize(1) && mask(pt(1),pt(2))==0)
        mask(pt(1),pt(2)) = i+1;
        qUsed = qUsed+1;
        q(:,qUsed) = [pt(1),pt(2), E(pt(1),pt(2))+pdist([pt(1),pt(2); py,px], distMethod)*distanceFac];
    end
    pt =[p(1)-1, p(2)];
    if(p(1) ~= 1 && mask(pt(1),pt(2))==0)
        mask(pt(1),pt(2)) = i+1;
        qUsed = qUsed+1;
        q(:,qUsed) = [pt(1),pt(2), E(pt(1),pt(2))+pdist([pt(1),pt(2); py,px], distMethod)*distanceFac];
    end
    pt =[p(1), p(2)+1];
    if(p(2) ~= Isize(2) && mask(pt(1),pt(2))==0)
        mask(pt(1),pt(2)) = i+1;
        qUsed = qUsed+1;
        q(:,qUsed) = [pt(1),pt(2), E(pt(1),pt(2))+pdist([pt(1),pt(2); py,px], distMethod)*distanceFac];
    end
    pt =[p(1), p(2)-1];
    if(p(2) ~= 1 && mask(pt(1),pt(2))==0)
        mask(pt(1),pt(2)) = i+1;
        qUsed = qUsed+1;
        q(:,qUsed) = [pt(1),pt(2), E(pt(1),pt(2))+pdist([pt(1),pt(2); py,px], distMethod)*distanceFac];
    end
   
    
    % doesnt look good when it draws around lower or right edge atm
    if(drawMode == 2 && mod(drawValue, 10)==0)
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

mask3d = repmat(mask,1,1,3);
Res(mask3d==i) = color(1,1,1);
mask3d(:,:,1)=0;
Res(mask3d==i) = color(1,1,2);
mask3d(:,:,2)=0;
Res(mask3d==i) = color(1,1,3);

if(drawMode == 1)
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


