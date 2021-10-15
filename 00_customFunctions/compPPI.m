% computes pixels/inch of a screen given its resolution (in pixels) & its
% diagonal (in inch)
function ppi = compPPI(widthPix, heightPix, diagInch)
    diagPix = sqrt(widthPix^2 + heightPix^2);
    ppi = diagPix/diagInch;
end