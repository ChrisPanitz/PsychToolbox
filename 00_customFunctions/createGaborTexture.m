function gaborTexture = createGaborTexture(window, sizePix, angle, spatFreqPix, baselineColor, peakColor)

sizePix = round(sizePix,0);

[x,y] = meshgrid(-sizePix(1)/2 : sizePix(1)/2, -sizePix(2)/2 : sizePix(2)/2); 
gaborFunction = (exp(-((x/sizePix(1)*4).^2)-((y/sizePix(2)*4).^2)) .* sin(cos(angle*pi/180)*(2*pi*spatFreqPix)*x + sin(angle*pi/180)*(2*pi*spatFreqPix)*y)); 
gaborPixels = baselineColor + peakColor * ((20./100)*gaborFunction);
% show it on screen 
gaborTexture = Screen('MakeTexture', window, gaborPixels, 0, 4);

end