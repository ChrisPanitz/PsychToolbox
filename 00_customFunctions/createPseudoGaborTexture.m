% creates a grating as texture object for Psychtoolbox on basis of a Gabor
% patch. Instead of smooth light-dark transitions, this grating's "bars" 
% have sharp edges. The grating has a circular shape.
% sizePix = [x, y] size in pixels
% angle = orientation of Grating in degrees with 0 == vertical orientation
% and clockwise rotation
% spatFreqPix = spatialFrequency of grating in units of cycles per pixel
% brightness = set of 3 values to determine brightness of stimulus components:
% [dark bars, background, bright bars]; e.g. [0 127 255] for black & white
% bars in front of gray background
function pseudoGaborTexture = createPseudoGaborTexture(window, sizePix, angle, spatFreqPix, brightness)

% to be safe, round size parameter
sizePix = round(sizePix,0);

% compute values for Gabor pattern; formula provided by Maeve Boylan (thanks!)
[x,y] = meshgrid(-sizePix(1)/2 : sizePix(1)/2, -sizePix(2)/2 : sizePix(2)/2); 
gaborFunction = (exp(-((x/sizePix(1)*4).^2)-((y/sizePix(2)*4).^2)) .* sin(cos(angle*pi/180)*(2*pi*spatFreqPix)*x + sin(angle*pi/180)*(2*pi*spatFreqPix)*y)); 

% "categorize" Gabor values into the dark, background and bright
gaborPixels = NaN(size(gaborFunction));
gaborPixels(gaborFunction < -.02) = brightness(1);
gaborPixels(gaborFunction >= -.02 & gaborFunction <= .02) = brightness(2);
gaborPixels(gaborFunction > .02) = brightness(3);

% create texture object
pseudoGaborTexture = Screen('MakeTexture', window, gaborPixels, 0, 4);
end