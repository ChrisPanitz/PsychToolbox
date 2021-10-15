function ppi = compPPI(widthPix, heightPix, diagInch)
    diagPix = sqrt(widthPix^2 + heightPix^2);
    ppi = diagPix/diagInch;
end