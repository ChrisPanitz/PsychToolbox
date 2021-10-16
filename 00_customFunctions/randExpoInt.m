% returns random interval duration in seconds, following an exponential
% function with mu as mean of the base function (exprnd function in MATLAB) 
% and a vector with minimal and maximum duration of interval; the value of
% 6.5 for computing lambda was obtained by trying out
function intInSec = randExpoInt(minMaxSec)
    intInSec = NaN;
    while ~(intInSec >= minMaxSec(1) && intInSec <= minMaxSec(2))
        lambda = 6.5/diff(minMaxSec);
        intInSec = minMaxSec(1) - log(1-rand(1,1))/lambda;
    end
end