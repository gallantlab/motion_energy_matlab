function params = preprocDownsample_GetMetaParams(Arg)
% Usage: params = preprocDownsample_GetMetaParams(Arg)
% 
% ML 2012.11.16

params.class = 'preprocDownsample';
switch Arg
    case 1
        % Simple box average, for TR=1, imhz=15
        params.dsType = 'box';
        params.imHz = 15;     % movie / image sequence frame rate
        params.sampleSec = 1; % TR
        params.frameshifts = []; % empty = no shift
        params.gaussParams = []; %[1,2]; % sigma,mean
    case 2
        % Simple box average, for TR=2, imhz=15
        params.dsType = 'box';
        params.imHz = 15;     % movie / image sequence frame rate
        params.sampleSec = 2; % TR
        params.frameshifts = []; % empty = no shift
        params.gaussParams = []; %[1,2]; % sigma,mean
    case 3
        % Gaussian downsampling, for TR=2, imhz=15
        params.dsType = 'gauss';
        params.imHz = 15;     % movie / image sequence frame rate
        params.sampleSec = 2; % TR
        params.frameshifts = []; % empty = no shift
        params.gaussParams = [1,2]; % mean, standard deviation
    case 4
        % Max downsampling, for TR=2, imhz=15
        params.dsType = 'max';
        params.imHz = 15;     % movie / image sequence frame rate
        params.sampleSec = 2; % TR
        params.frameshifts = []; % empty = no shift
        params.gaussParams = []; %[1,2]; % sigma,mean
    case 5
        % Simple box average, for TR=2, imhz=24
        params.dsType = 'box';
        params.imHz = 24;     % movie / image sequence frame rate
        params.sampleSec = 2; % TR
        params.frameshifts = []; % empty = no shift
        params.gaussParams = []; %[1,2]; % sigma,mean
    case 6
        % Simple box average, for TR=1, imhz=24
        params.dsType = 'box';
        params.imHz = 24;     % movie / image sequence frame rate
        params.sampleSec = 1; % TR
        params.frameshifts = []; % empty = no shift
        params.gaussParams = []; %[1,2]; % sigma,mean        
    otherwise
        error('Unknown parameter configuration!');
end
