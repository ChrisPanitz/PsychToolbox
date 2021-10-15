% wait for any mouse click, detected via the GetMouse function in
% Psychtoolbox
function [x, y, buttons] = waitForClick
    buttons = 0;
    while ~any(buttons) % wait for press
        [x, y, buttons] = GetMouse;
    
        % Wait 10 ms before checking the mouse again to prevent
        % overload of the machine at elevated Priority()
        WaitSecs(0.01);
    end
end