function intInSec = randExpoInt(mu, minMaxSec)
    intInSec = NaN;
    
    while ~(intInSec >= minMaxSec(1) && intInSec <= minMaxSec(2))
        intInSec = minMaxSec(1) + exprnd(mu,1,1)*diff(minMaxSec);
    end
end
