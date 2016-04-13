function params = preprocDownsample_GetMetaParams(Arg)
% Usage: params = preprocDownsample_GetMetaParams(Arg)
% 
% ML 2012.11.16

params.class = 'preprocDownsample';
switch Arg
    case 1
        % Simple box average
        params.dsType = 'box';
        params.imHz = 15;     
        params.sampleSec = 2; 
        params.frameshifts = []; % empty = no shift
        params.gaussParams = []; %[1,2]; % sigma,mean
    case 2
        % Gaussian downsampling
        params.dsType = 'gauss';
        params.imHz = 15;
        params.sampleSec = 2;
        params.frameshifts = []; % empty = no shift
        params.gaussParams = [1,2]; % mean, standard deviation
    case 3
        % Max downsampling 
        params.dsType = 'max';
        params.imHz = 15;
        params.sampleSec = 2;
        params.frameshifts = []; % empty = no shift
        params.gaussParams = []; %[1,2]; % sigma,mean

    case 4
        % Simple box average, w/ different stimulus presentation rate
        params.dsType = 'box';
        params.imHz = 24;     
        params.sampleSec = 2; 
        params.frameshifts = []; % empty = no shift
        params.gaussParams = []; %[1,2]; % sigma,mean
    case 101
        % SHADY: imHz SHOULD be read in from stimulus input, not specified
        % as a parameter here. This is just a kloodge for the HCP data,
        % for which the movies were presented at a different rate than we
        % usually present (24 hz)
        % Simple box average
        params.dsType = 'box';
        params.imHz = 24;     % These two values will be overwritten by stimulus 
        params.sampleSec = 1; % parameters in Stimulus class input to preprocDownsample
        params.frameshifts = []; % empty = no shift
        params.gaussParams = []; %[1,2]; % sigma,mean        
    otherwise
        error('Unknown parameter configuration!');
end
