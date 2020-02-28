clc
clear all

% Find boundaries for CIELab based on representable RGB

color = [];
for r = 0:0.1:1
    for g = 0:0.1:1
        for b = 0:0.1:1
            color = [color; r, g, b];            
        end
    end
end

lab = rgb2lab(color);
maxLab = max(lab);
minLab = min(lab);