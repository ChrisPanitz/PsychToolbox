% returns random interval duration in seconds, following an exponential
% function with mu as mean of the base function (exprnd function in MATLAB) 
% and a vector with minimal and maximum duration of interval; mu = .1 works well
function intInSec = randExpoInt(mu, minMaxSec)
    intInSec = NaN;
    while ~(intInSec >= minMaxSec(1) && intInSec <= minMaxSec(2))
        intInSec = minMaxSec(1) + exprnd(mu,1,1)*diff(minMaxSec);
    end
end