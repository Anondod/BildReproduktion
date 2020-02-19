% expansion test

% 0 draws nothing
% 1 draws every 1% of total segments finished
% 2 draws when each segment is finished
% 3 shows every pixel update
drawMode = 1;

% have 1px contours around segments?
doContours = 1;

% for step 2 and the total coverage at end
drawcol(1,1,:) = [0,0,1];

distCityBlock = @(p1,p2) sum(abs(p2-p1));
distEuclidean = @(p1,p2) sqrt(sum((p2-p1).^2));

% choose one from above: euclidean , cityblock
distFun = distCityBlock;
distanceFac = 0.4;

qSize = 10000;

maxErr = 10000;

% 0 - choose at random.
% anything else - random place of the first
% firstSelection black mask positions.
firstSelection = 100;

% if a pixel of the segment has to large an error the segment stops
cancelThreshold = 80;

%nr mosaics (max)
n = 5000;

% Max Segment radius size
SEGMaxRad = 50;


% DECIDE INPUT IMAGE HERE
I = im2double(imresize(imread('test3.jpg'),1));

% Blur image?
%I = imgaussfilt(I,2);
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

% alt faster ver. but looks worse maybe. also breaks coverage.
% firstSelection = 1 gives the same result every time (maybe useful)
if(firstSelection ~= 0)
    [py_temp,px_temp] = find(mask==0,1,'first');
else
    [py_temp,px_temp] = find(mask==0);
end

unfilled = size(px_temp,1);

if(unfilled==0)
    break;
end

index = randi(unfilled);
py = py_temp(index);
px = px_temp(index);


%progress bar
if(mod(i,printmod)==0)
    if(firstSelection~=0)
        covered = (px*Isize(1)+py) /(Isize(1)*Isize(2));
    else
        covered = 1-unfilled/(Isize(1)*Isize(2));
    end
    fprintf(i/2+"/"+n+" ("+100*i/(2*n)+"%%) segments done\n");
    fprintf(100*(covered)+"%% coverage\n\n");
    
    if(drawMode == 1)
        warning off;
        imshow(Res)
        warning on;
    end
    drawnow
end


color = I(py,px,:);

SEGxstart = max(1,px-SEGMaxRad);
SEGxend = min(Isize(2),px+SEGMaxRad);
SEGystart = max(1,py-SEGMaxRad);
SEGyend = min(Isize(1),py+SEGMaxRad);

SEGx = SEGxstart:SEGxend;
SEGy = SEGystart:SEGyend;

SEGsize(1) = size(SEGy,2);
SEGsize(2) = size(SEGx,2);

e = Ilab(SEGy,SEGx,:) - Ilab(py,px,:);
E = sqrt(e(:,:,1).^2+e(:,:,2).^2+e(:,:,3).^2);

SEGmask = mask(SEGy,SEGx);

px = px - (max(1,px-SEGMaxRad)-1);
py = py - (max(1,py-SEGMaxRad)-1);

qUsed = 1;
q(:,1) = [py px 0];

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
        
    SEGmask(p(1),p(2)) = i;
    
    % adding adjacents to queue and mask
    
    % pdist substitute: (sqrt((px-pt(2)).^2+(py-pt(1)).^2))

    
    pt =[p(1)+1, p(2)];
    if(p(1) ~= SEGsize(1) && SEGmask(pt(1),pt(2))==0)
        SEGmask(pt(1),pt(2)) = i+1;
        qUsed = qUsed+1;
        q(:,qUsed) = [pt(1),pt(2), E(pt(1),pt(2))+distFun([pt(1),pt(2)],[py,px])*distanceFac];
    end
    pt =[p(1)-1, p(2)];
    if(p(1) ~= 1 && SEGmask(pt(1),pt(2))==0)
        SEGmask(pt(1),pt(2)) = i+1;
        qUsed = qUsed+1;
        q(:,qUsed) = [pt(1),pt(2), E(pt(1),pt(2))+distFun([pt(1),pt(2)],[py,px])*distanceFac];
    end
    pt =[p(1), p(2)+1];
    if(p(2) ~= SEGsize(2) && SEGmask(pt(1),pt(2))==0)
        SEGmask(pt(1),pt(2)) = i+1;
        qUsed = qUsed+1;
        q(:,qUsed) = [pt(1),pt(2), E(pt(1),pt(2))+distFun([pt(1),pt(2)],[py,px])*distanceFac];
    end
    pt =[p(1), p(2)-1];
    if(p(2) ~= 1 && SEGmask(pt(1),pt(2))==0)
        SEGmask(pt(1),pt(2)) = i+1;
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
        
        draw = I(scope1,scope2,:) + double(SEGmask(scope1,scope2) == i).*drawcol;
        %draw = E(scope1,scope2,:)/100 + double(mask(scope1,scope2) == i).*drawcol;
        imshow(draw);
        
        %q(:,1:qUsed+1) % what's in q at the moment
        pause(0.001);
    end
    
    drawValue = drawValue+1;
end

if(doContours==0)
    SEGmask(SEGmask==i+1) = 0;
end

mask(SEGy,SEGx) = SEGmask;

SEGmask3d = repmat(SEGmask,1,1,3);

TempIM = Res(SEGy,SEGx,:);
TempIM(SEGmask3d==i) = color(1,1,1);
SEGmask3d(:,:,1)=0;
TempIM(SEGmask3d==i) = color(1,1,2);
SEGmask3d(:,:,2)=0;
TempIM(SEGmask3d==i) = color(1,1,3);

Res(SEGy,SEGx,:) = TempIM;

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
