function [flickVec] = genFlickVec(frPerCyc, frTotal, func)
    if strcmp(func, 'box')
        flickVec = repmat([ones(frPerCyc/2,1); zeros(frPerCyc/2,1)], ...
                          frTotal/frPerCyc, 1)';
    elseif strcmp(func, 'sin')
        flickVec = repmat(sin(1.5*pi:(2*pi/frPerCyc):(3.5*pi-2*pi/frPerCyc)) / 2 + .5, ...
                          1, frTotal/frPerCyc);
    end
end


