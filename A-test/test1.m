clear all;
I = im2double(imread('elin.jpg'));

th = im2double(imread('thread5.png'));
xth = size(th,2);
yth = size(th,1);
%xstart = ceil(xth/2);
%ystart = ceil(yth/2);


xseg = 100;


yseg = floor(xseg * size(I,1)/size(I,2) * xth/yth);
xstep = size(I,2)/xseg;
ystep = size(I,1)/yseg;

stitchcolors = zeros(yseg,xseg,3);
for i=1:xseg
    i1 = floor((i-1)*xstep+1);
    i2 = floor(i1+xstep-1);
    for j=1:yseg
        j1 = floor((j-1)*ystep+1);
        j2 = floor(j1+ystep-1);
        stitchcolors(j,i,:) = mean(mean(I(j1:j2,i1:i2,:)));
    end
end


imshow(stitchcolors)
pause




clear a
a = ones(yth*yseg,xth*xseg,3) * 255;
%a(ystart:yspace:end,xstart:xspace:end,:) = 100;

maskfac = 1;
thmask = repmat((mean(th,3)<=maskfac), 1, 1, 3);

for i=1:xseg
    i2 = (i-1)*xth+1;
    for j=1:yseg
        j2 = (j-1)*yth+1;
        temp = th*0.8;
        temp_c = th.*stitchcolors(j,i,:);
        temp(thmask)=temp_c(thmask);
        a(j2:j2+yth-1,i2:i2+xth-1,:) = temp;
    end
end

imshow(a);
imwrite(a,'a.png')
