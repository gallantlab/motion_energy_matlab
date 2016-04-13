function params = preprocNormalize_GetMetaParams(argNum)
% Usage: params = preprocNormalize_GetMetaParams(argNum)
% 
% Get meta-parameters for preprocNormalize
% 
% ML 2013.03.21

params.class = 'preprocNormalize';
switch argNum
    case 1
        % Original arguments recmomended by SN
        params.valid_w_index = []; % specific index of channels to keep (overrides .reduceChannels)
        params.reduceChannels = []; % number of channels / pct of channels to keep
        params.normalize = 'zscore'; % normalization method
        params.crop = []; % min/max to which to crop; empty does nothing
    case 2
        params.valid_w_index = [];
        params.reduceChannels = []; %
        params.normalize = 'gaussianize';
        params.crop = [];
    case 3
        % Original arguments recmomended by SN; crops to [-3.5,3.5]
        params.valid_w_index = []; % specific index of channels to keep (overrides .reduceChannels)
        params.reduceChannels = []; % number of channels / pct of channels to keep
        params.normalize = 'zscore'; % normalization method
        params.crop = [-3.5,3.5]; % min/max to which to crop; empty does nothing
    case 4
        % blank for now...
    case 5
        % Original arguments recmomended by SN; crops to [-5,5]
        params.valid_w_index = []; % specific index of channels to keep (overrides .reduceChannels)
        params.reduceChannels = []; % number of channels / pct of channels to keep
        params.normalize = 'zscore'; % normalization method
        params.crop = [-5,5]; % min/max to which to crop; empty does nothing
        params.useTrnParams = true;
end