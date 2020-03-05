function res = centroid2PaletteColor(centroid,selectionColor)
% finds the nearest palette color from selectionColor for each centroid
% (with no duplicate results)

selcol = selectionColor;

res = centroid;
for i = 1:length(centroid)
    a = sqrt(sum((selcol-centroid(i,:)).^2,2));
    index = find(min(a) == a , 1, 'first');
    res(i,:) = selcol(index,:);
    selcol(index,:) = [];
end
    
end

