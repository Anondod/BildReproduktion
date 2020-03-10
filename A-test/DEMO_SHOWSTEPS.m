warning off;
clear all;
warning on;


% **INTERESTING OPTIONS SECTION START**

    % DECIDE INPUT IMAGE HERE
    I = im2double(imresize(imread('Images/elin.jpg'),1));

    % Blur image?
    %I = imgaussfilt(I,2);
    

    % 0 draws nothing
    % 1 draws every 1% of total segments finished
    % 2 draws when each segment is finished
    % 3 shows every pixel update
    drawMode = 3;

    % 0 - choose at random.
    % 1 - choose first open spot
    % 2 - choose most monotone areas first (low derivative)
    selectionMethod = 0;
    
    % 0 - all colors
    % 1 - colors from palette
    % 2 - image optimal colors from palette
    colorMode = 0;
    
    % if false - segments are created normally then colored
    % if true - segments are depending on the palette color
    compareWithPalette = true;
    
    palette = load('palette6.mat');
    palette = palette.a;
    % Number of colors to take from palette with colorMode 2
    n_clusteringColors = 20;
    

    % Use morphology on the segment + contours(if true)
    doMorph = true;
    % Close structuring element (STREL) 
    closeMorph1 = strel('line', 20, 90);
    closeMorph2 = strel('line', 20, 0);
    closeMorph3 = strel('sphere', 4);

    doContours = false;
    contoursMorph = strel('sphere', 3);

    
    distCityBlock = @(p1,p2) sum(abs(p2-p1));
    distEuclidean = @(p1,p2) sqrt(sum((p2-p1).^2));
    % choose one from above: euclidean , cityblock
    distFun = distEuclidean;
    distanceFactor = 0.3;

    % allowed accumulated error per segment
    maxErr = 1000;
    % if a pixel of the segment has to large an error the segment stops
    cancelThreshold = 20;
    % Max Segment section radius
    SEGMaxRad = 40;

    
    %nr mosaics (max)
    n = 5000;
    % max segment size(pixels)
    qSize = 10000;

% **INTERESTING OPTIONS SECTION END**


Isize = [size(I,1) size(I,2)];
Ilab = rgb2lab(I);

% Dmag used for selectionMode 2 (image derivative)
[magR,dir] = imgradient(I(:,:,1)); 
[magG,dir] = imgradient(I(:,:,2)); 
[magB,dir] = imgradient(I(:,:,3));
Dmag = sqrt(magR.^2 + magG.^2 + magB.^2);
Dmag = imgaussfilt(Dmag/max(max(Dmag)),4);

% initialization of data structures
Res = zeros(size(I));
mask = uint16(zeros(Isize(1),Isize(2)));
q = zeros(4, qSize);

% how man segments between "progress bar" prints
printmod = ceil(n/50);
printmod = printmod + mod(printmod,2);

% for step 2 and the total coverage at end
drawcol(1,1,:) = [0,0,1];
drawcol2(1,1,:) = [1,0,0];

covered = 0;



if(colorMode == 2)
    disp("Starting clustering...");
    drawnow;
    
    L = Ilab(:,:,1);
    a = Ilab(:,:,2);
    b = Ilab(:,:,3);

    img(:,1) = L(:);
    img(:,2) = a(:);
    img(:,3) = b(:);

    [clusterInd, theLABColors] = kmeans(img, n_clusteringColors);
    
    disp("Done clustering");
    drawnow;
    
    palette = rgb2lab(palette);
    palette = centroid2PaletteColor(theLABColors, palette);
    
elseif(colorMode==1)
    palette = rgb2lab(palette);
end


% main segment for-loop
for i=2:2:2*n

% alt faster ver. but looks worse maybe. also breaks coverage.
% firstSelection = 1 gives the same result every time (maybe useful)
if(selectionMethod == 0)
    [py_temp,px_temp] = find(mask==0);
    
    unfilled = size(px_temp,1);
    if(unfilled==0)
        break;
    end
    
    index = randi(unfilled);
    py = py_temp(index);
    px = px_temp(index);
    
elseif(selectionMethod == 1)
    [py,px] = find(mask==0,1,'first');
    
    if(size(px,1)==0)
        break;
    end
    
else %selectionMethod == 2
    [py,px] = find(Dmag==min(min(Dmag)));
    
    [py_temp,px_temp] = find(mask==0);
    unfilled = size(px_temp,1);
        
    if(unfilled==0)
        break;
    end
    
    px = px(1);
    py = py(1);
end

if(colorMode~=0)
    % color selection
    for j = 1:length(palette)
        test(:) = Ilab(py,px,:);
        e = palette(j, :) - test;
        ETheColors(j, :) = sqrt(e(1).^2+e(2).^2+e(3).^2);
    end
    LABcolor = palette(ETheColors == min(ETheColors), :);
    color(1,1,:) = lab2rgb(LABcolor)';
else
    color = I(py,px,:);
end

%progress bar
if(mod(i,printmod)==0)
    if(selectionMethod==0)
        covered = 1-unfilled/(Isize(1)*Isize(2));
    elseif(selectionMethod==1)
        covered = (px*Isize(1)+py) /(Isize(1)*Isize(2));
    else %selectionMethod == 2
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


SEGxstart = max(1,px-SEGMaxRad);
SEGxend = min(Isize(2),px+SEGMaxRad);
SEGystart = max(1,py-SEGMaxRad);
SEGyend = min(Isize(1),py+SEGMaxRad);

SEGx = SEGxstart:SEGxend;
SEGy = SEGystart:SEGyend;

SEGsize(1) = size(SEGy,2);
SEGsize(2) = size(SEGx,2);

if(colorMode~=0 && compareWithPalette)
    comparisonLabColor(1,1,:) = LABcolor';
else
    comparisonLabColor = Ilab(py,px,:);
end
e = Ilab(SEGy,SEGx,:) - comparisonLabColor;
E = sqrt(e(:,:,1).^2+e(:,:,2).^2+e(:,:,3).^2);

SEGmask = mask(SEGy,SEGx);

px = px - (max(1,px-SEGMaxRad)-1);
py = py - (max(1,py-SEGMaxRad)-1);

qUsed = 1;
q(:,1) = [py px 0 0];

drawValue = 0;
% Loop per segment
acc_error = 0;
while acc_error < maxErr
    [err,lowest] = min(q(3,1:qUsed));
    p = q(1:2,lowest);
    acc_error = acc_error + err;
    
    colErr = q(4,lowest);
    
    % break if q is empty or
    if(qUsed ~= 0)
        q(:,lowest) = q(:,qUsed);
        qUsed = qUsed-1;
        if(colErr > cancelThreshold)
            continue;
        end
    else
        break;
    end
        
    SEGmask(p(1),p(2)) = i;
    
    % ADD ADJACENT PIXELS TO q AND mask
    pt =[p(1)+1, p(2)];
    if(p(1) ~= SEGsize(1) && SEGmask(pt(1),pt(2))==0)
        SEGmask(pt(1),pt(2)) = i+1;
        qUsed = qUsed+1;
        q(:,qUsed) = [pt(1),pt(2), E(pt(1),pt(2))+distFun([pt(1),pt(2)],[py,px])*distanceFactor, E(pt(1),pt(2))];
    end
    pt =[p(1)-1, p(2)];
    if(p(1) ~= 1 && SEGmask(pt(1),pt(2))==0)
        SEGmask(pt(1),pt(2)) = i+1;
        qUsed = qUsed+1;
        q(:,qUsed) = [pt(1),pt(2), E(pt(1),pt(2))+distFun([pt(1),pt(2)],[py,px])*distanceFactor, E(pt(1),pt(2))];
    end
    pt =[p(1), p(2)+1];
    if(p(2) ~= SEGsize(2) && SEGmask(pt(1),pt(2))==0)
        SEGmask(pt(1),pt(2)) = i+1;
        qUsed = qUsed+1;
        q(:,qUsed) = [pt(1),pt(2), E(pt(1),pt(2))+distFun([pt(1),pt(2)],[py,px])*distanceFactor, E(pt(1),pt(2))];
    end
    pt =[p(1), p(2)-1];
    if(p(2) ~= 1 && SEGmask(pt(1),pt(2))==0)
        SEGmask(pt(1),pt(2)) = i+1;
        qUsed = qUsed+1;
        q(:,qUsed) = [pt(1),pt(2), E(pt(1),pt(2))+distFun([pt(1),pt(2)],[py,px])*distanceFactor, E(pt(1),pt(2))];
    end
   
    
    % draw per pixel (or fewer depending on mod value)
    if(drawMode == 3 && mod(drawValue, 1)==0)
        draw = I(SEGy,SEGx,:);
        channel = repmat(SEGmask, 1,1,3);
        channel(:,:,2:3) = 0;
        draw(channel == i) = 1;
        channel(:,:,1) = 0;
        channel(:,:,2) = SEGmask;
        draw(channel == (i+1)) = draw(channel == (i+1)) + 0.4;
        imshow(imresize(draw, 10, 'nearest'));
        pause(0.01);
    end
    
    drawValue = drawValue+1;
end


if(doContours==false)
    SEGmask(SEGmask==i+1) = 0;
end

draw = I(SEGy,SEGx,:);
channel = repmat(SEGmask, 1,1,3);
channel(:,:,2:3) = 0;
draw(channel == i) = 1;
channel(:,:,1) = 0;
channel(:,:,2) = SEGmask;
draw(channel == (i+1)) = 1;
imshow(imresize(draw, 10, 'nearest'));

if(doMorph == true)
    tempClose = SEGmask == i;

    tempClose = imclose(tempClose, closeMorph1);
    tempClose = imclose(tempClose, closeMorph2);
    tempClose = imclose(tempClose, closeMorph3);


    if(doContours==true)
        tempDilate = imdilate(tempClose, contoursMorph);
        contorsResult = logical(tempDilate - tempClose);
        SEGmask(contorsResult) = i+1;
    end

    SEGmask(tempClose) = i;
end

if(selectionMethod==2)
    Dmag(SEGy,SEGx) = Dmag(SEGy,SEGx)+double(SEGmask);
end

pause;
draw = I(SEGy,SEGx,:);
channel = repmat(SEGmask, 1,1,3);
channel(:,:,2:3) = 0;
draw(channel == i) = 1;
channel(:,:,1) = 0;
channel(:,:,2) = SEGmask;
draw(channel == (i+1)) = 1;
imshow(imresize(draw, 10, 'nearest'));
pause;

mask(SEGy,SEGx) = SEGmask;

SEGmask3d = repmat(SEGmask,1,1,3);

TempIM = Res(SEGy,SEGx,:);
TempIM(SEGmask3d==i) = color(1,1,1);
SEGmask3d(:,:,1)=0;
TempIM(SEGmask3d==i) = color(1,1,2);
SEGmask3d(:,:,2)=0;
TempIM(SEGmask3d==i) = color(1,1,3);

SEGmask3d = repmat(SEGmask,1,1,3);

if(doContours == true)
    TempIM(SEGmask3d==i+1) = 0;
    mask3d(:,:,1)=0;
    TempIM(SEGmask3d==i+1) = 0;
    mask3d(:,:,2)=0;
    TempIM(SEGmask3d==i+1) = 0;
end

Res(SEGy,SEGx,:) = TempIM;

if(drawMode == 2)
    warning off;
    imshow(Res)
    warning on;
end

end

fprintf("done at " + 100*(1-unfilled/(Isize(1)*Isize(2)))+"%% coverage\n");

%%
warning off;

close all
figure;
draw = I;
imshow(draw);
title("original");

figure;
imshow(double(mask)/i)
title("segments");


figure;
imshow(Res)
title("result");

warning on;

% SAVE COMMAND
%imwrite(Res,'insertname.png')
