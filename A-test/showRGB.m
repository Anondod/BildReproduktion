function showRGB(RGB)
% showRGB(RGB)
%
% RGB = Matris med RGB-värden som ska visas på skärmen.
% Matrisen bör se ut som följande:
% RGB=[r1 g1 b1;r2 g2 b2;...]

% Martin Solli, marso@itn.liu.se
% ITN, Linköpings Universitet

r=size(RGB,1);

colors=[];
co=0;
ro=0;
for n=1:r
    if (co==10) co=0; ro=ro+1; end;
    colors(ro*75+1:ro*75+75,co*75+1:co*75+75,1)=RGB(n,1);
    colors(ro*75+1:ro*75+75,co*75+1:co*75+75,2)=RGB(n,2);
    colors(ro*75+1:ro*75+75,co*75+1:co*75+75,3)=RGB(n,3);
    co=co+1;
end;

figure;
imshow(colors);