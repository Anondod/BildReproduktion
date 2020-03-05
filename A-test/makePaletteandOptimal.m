function [opt_colors, palette] = makePaletteandOptimal(steps, n_opt_colors)
% returns a palette of steps^3 colors (~total pearls) and an 'optimized'
% palette from this palette

step = 1/(steps-1);

palette = [];
for r = 0:step:1
    for g = 0:step:1
        for b = 0:step:1
            palette = [palette; r, g, b];            
        end
    end
end

labColors = rgb2lab(palette);

[~, optLab] = kmeans(labColors, n_opt_colors-2);

optLab = centroid2PaletteColor(optLab, labColors);

% add black and white (almost never included from the clustering)
optLab = [optLab; 0, 0, 0];
optLab = [optLab; 100, 0, 0];

opt_colors = lab2rgb(optLab);
opt_colors(opt_colors<=0.00000000001) = 0;
opt_colors(opt_colors>1) = 1;

showRGB(palette)
showRGB(opt_colors)
end

