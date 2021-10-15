function nrPix = pixFromAngle(visAngle, distInch, ppi)
    % visAngle = 2 * atan(imSize/(2dist)) * 180/pi
    % visAngle = 2 * atan(nrPix/ppi/(2dist)) * 180/pi
    % visAngle/2/180*pi = atan(nrPix/ppi/2/dist)
    % tan(visAngle/360*pi) = nrPix/ppi/2/dist
    % tan(visAngle/360*pi) = nrPix/(ppi*2*dist)
    % nrPix = tan(visAngle/360*pi) * ppi * 2 * dist
   
    % that's actually the image size in units
    % nrPix = tan(visAngle/360*pi) * 2 * distInch;
    % that's actual pixels
    nrPix = round(tan(visAngle/360*pi) * ppi * 2 * distInch, 0);
end

    % visAngle = 2 * atan(imSize/(2dist))
    % visAngle/2 = atan(imSize/(2dist))
    % tan(visAngle/2) = imSize/(2dist)
    % imSize = tan(visAngle/2) * 2 * dist
