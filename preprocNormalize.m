function varargout = preprocNormalize(S,params)
% Usage: varargout = preprocNormalize(S,params)
% 
% Normalize channels of a model. Most commonly, z-score each channel, but
% other methods are available.
% 
% Inputs: 
%   S : stimulus / preprocessed stimulus. Better be 2D (time x channels)
%   params : parameter struct array, with fields: 
%       .normalize : string, one of the following: 'zscore' [default],
%               'gaussianize', 'uniform', '0to1', '-1to1' 
%       .reduceChannels : if scalar < 1, keep all channels with stds >
%               params.reduceChannels * max std; if scalar > 1, keep n
%               channels; if 1 or true, use following parameter as index to
%               keep some channels. Default = [] = do nothing.
%       [.reduceChannelsValidChannels] : index of channels to keep
%               (optional, only use if .reduceChannels==1)
%       .crop : 2-element vector [min,max] - crop values above/below
%               max/min to max/min. Default = [] = do nothing.
% 
% Outputs: 
%   Spreproc : normalized stimulus
%   [params] : preprocessing params, potentially w/ .means, .stds added (if
%           params.normalize = 'zscore')
% 
% ML 2013.03.20

% TO DO: incorporate sparseness measures per channel? Is this even useful?

% Default parameters
dParams.class = 'preprocNormalize';
dParams.reduceChannels = []; 
dParams.crop = []; 
dParams.normalize = 'zscore'; % 'gaussianize','makeUniform','0to1','-1to1' (?)
% Fill params w/ defaults
if ~exist('params','var')
    params = struct;
end
params = defaultOpt(params,dParams);
% Return params if no input
if ~nargin
    varargout{1} = dParams;
    return
end

% Normalize
switch params.normalize
    case 'zscore'
        if isfield(params, 'means') % Already preprocessed. Use the means and stds
            [Spreproc] = norm_std_mean(S, params.stds, params.means);
        else
            [Spreproc, stds, means] = norm_std_mean(S);
            params.means = means;
            params.stds = stds;
        end
    case 'gaussianize'
        Spreproc = zeros(size(S));
        for iCh = 1:size(S,2)
            SppU = makeUniform(S(:,iCh));
            Spreproc(:,iCh) = norminv(SppU);
        end
    case 'uniform'
        Spreproc = zeros(size(S));
        for iCh = 1:size(S,2)
            Spreproc(:,iCh) = makeUniform(S(:,iCh));
        end
    case '0to1'
        Sr = bsxfun(@minus,S,min(S,[],1));
        Spreproc = bsxfun(@rdivide,Sr,max(abs(Sr),[],1));
    case '-1to1'
        Spreproc = bsxfun(@rdivide,S,max(abs(S),[],1));        
end

% Reduce number of channels based on the standard deviation of channels, or
% a pre-defined index (reduceChannelsValidChannels) 
if ~isempty(params.reduceChannels)
    warning('Reducing channels is not well tested in re-implementation of code yet!')
    orig_nch = size(Spreproc,2);
    if isfield(params, 'reduceChannelsValidChannels') % Already preprocessed. Use valid channels.
        v = params.reduceChannelsValidChannels;
        Spreproc = Spreproc(:,v);
    else
        if ~isfield(params, 'stds')
            [d, stds, means] = norm_std_mean(Spreproc);
        else
            stds = params.stds;
        end
        if params.reduceChannels < 1
            maxstd = max(stds);
            v = find(stds >= maxstd*params.reduceChannels);
            Spreproc = Spreproc(:,v);
        else
            [d s] = sort(stds, 'descend');
            v = s(1:min([length(s) params.reduceChannels]));
            v = sort(v);
            Spreproc = Spreproc(:,v);
        end
        params.reduceChannelsValidChannels = v;
    end
    if isfield(params,'verbose') && params.verbose, fprintf('Trucated channels %d -> %d.\n', orig_nch, length(v)); end
end

% Crop values to fit within specified range
if ~isempty(params.crop)
    Spreproc = max(Spreproc,params.crop(1));
    Spreproc = min(Spreproc,params.crop(2));
end
% Output
varargout{1} = Spreproc;
if nargout>1
    varargout{2} = params;
end
