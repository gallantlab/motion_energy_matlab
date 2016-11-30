function varargout = preprocDownsample(S, params)
% Usage: [Spp, params] = preprocDownsample(S, params)
% 
% Downsamples stimuli to (fMRI or other) sampling rate; e.g., takes frames
% of a movie presented at 15 Hz and downsamples them to an fMRI sampling
% rate (TR) of .5 Hz (2 seconds / measurement). 
% 
% Inputs: 
%   S = preprocessed stimulus matrix, [time x channels]
%   params = parameter struct, with fields:
%       .dsType = string specifying type of downsampling: 'box' [default],
%           'gauss', 'max', or 'none'
%       .gaussParams = 2-element array specifying [,] (??). Only necessary
%           if params.dsType = 'gauss'
%       .imHz = frame rate of stimulus in Hz 
%       .sampleSec = length in seconds of 1 sample of data (for fMRI, this
%           is the TR or repetition time)
%       .frameshifts = amount to shift frames; empty implies no shift
%       .gaussParams = standard deviation and temporal offset for Gaussian
%           downsampling window
% Output:
%   Spp

% Default parameters
dParams.dsType = 'box';
dParams.imHz = 15;
dParams.sampleSec = 2;
dParams.frameshifts = []; % empty = no shift
dParams.gaussParams = []; %[1,2]; % sigma, offset
% Fill in default params
if ~exist('params','var')
    params = struct;
end
params = defaultOpt(params,dParams);
% Return params if no inputs
if ~nargin
    varargout{1} = params;
    return
end

fr_per_sample = params.sampleSec*params.imHz;
% downsample the preprocessed stimuli
switch params.dsType
    case 'box'
        if isfield(params,'frameshifts') && ~isempty(params.frameshifts)
          fprintf('shifting %d frames...\n', params.frameshifts);
          S=circshift(S,[params.frameshifts 0]);
        end
        tframes = floor(size(S,1)/fr_per_sample)*fr_per_sample;
        S = S(1:tframes,:);
        S = reshape(S, fr_per_sample, [], size(S,2));
        S = reshape(mean(S,1), [], size(S,3));
    case 'none'
        0; % do nothing
    case 'max'
        tframes = floor(size(S,1)/fr_per_sample)*fr_per_sample;
        S = S(1:tframes,:);
        S = reshape(S, fr_per_sample, [], size(S,2));
        S = reshape(max(S,[],1), [], size(S,3));
    case 'gauss'
        ksigma = params.gaussParams(1);
        if ksigma~=0
            ki = -ksigma*2.5:1/fr_per_sample:ksigma*2.5;
            k = exp(-ki.^2/(2*ksigma^2));
            S = conv2(S, k'/sum(k), 'same');
        end
        sonset = 7;
        if length(params.gaussParams)>=2
            sonset = params.gaussParams(2);
        end
        S = S(sonset:fr_per_sample:end,:);
end

% Output
varargout{1} = S;
if nargout>1
    varargout{2} = params;
end
