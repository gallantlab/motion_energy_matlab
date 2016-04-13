function params = preprocNonLinearOut_GetMetaParams(argNum)
% Usage: params = preprocNonLinearOut_GetMetaParams(argNum)
% 
% Get preset arguments for preprocNonLinearOut
%
% 

params.class = 'preprocNonLinearOut';
switch argNum
    case 1
        % Original params recommended by SN
        params.gainControl = []; % Broken: gain control for each channel based on luminance /color.
        params.gainControlOut = []; % Broken: gain control for each channel based on luminance /color.
        params.nonLinOutExp = 'log'; % Output nonlinearity
        params.nonLinOutParam = 1.0000e-05; % delta to add to channel values to prevent log(0) = -inf
    case 2
        % Original params recommended by SN
        params.gainControl = []; % Broken: gain control for each channel based on luminance /color.
        params.gainControlOut = []; % Broken: gain control for each channel based on luminance /color.
        params.nonLinOutExp = .5; % Output nonlinearity
        %params.nonLinOutParam = 1.0000e-05; % delta to add to channel values to prevent log(0) = -inf
        
end
