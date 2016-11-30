function varargout = preprocNonLinearOut(S, params)
% Usage: [Spreproc,params] = preprocNonLinearOut(S, params)
% 
% Output nonlinearity on each channel of a model. Generally done BEFORE
% (zscore or other) normalization. 
% 
% Inputs: 
%   S : Stimulus or PreprocessedStimulus (STRFlab classes), or simple
%     numerical matrix (assumed to be nFrames x nChannels)
%   params : a struct array of parameters, with fields: 
%       .nonLinOutExp : exponent to which to raise each channel value,
%                       thus: abs(S).^x .*sign(S)  Multiple values (e.g.
%                       [.5,2]) raise each column to each different
%                       exponent, and concatenate the results (here,
%                       doubling the number of channels). 
%                       Alternately, this can be 'log', w/ another
%                       parameter (see next)
%       .nonLinOutParam = []; % This needs filling in for certain options
%                       of nonLinOutExp (e.g. for 'log', it specifies a
%                       small delta to add to avoid log(0)= -inf
%       .verbose = T/F, verbose printing
% 
% Modified from SN code by ML 2013.03.21

% Default parameters
dParams.class = 'preprocNonLinearOut';
dParams.nonLinOutExp = 0.5; % Exponent to which to raise each value
dParams.verbose = true;
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

if isfield(params,'verbose') && params.verbose
    fprintf('Processing static output nonlinearity...'); 
end
% .nonLinOutParam has a different meaning for each different output nonlinearity.
if ischar(params.nonLinOutExp)
    switch params.nonLinOutExp
        case {'log'}
            % For log, d is a small value to add to assure no -Inf channels.
            d = params.nonLinOutParam;
            for ii=1:size(S,2)
                S(:,ii) = log(d+S(:,ii));
            end
        case {'logstd'}
            if isfield(params,'nonLinOutStd')
                stds = params.nonLinOutStd;
            else
                stds = nanstd(S);
                params.nonLinOutStd = stds;
            end
            % For logstd, d is also a small value to avoid -Inf
            d = params.nonLinOutParam;
            for ii=1:size(S,2)
                S(:,ii) = log(d+S(:,ii)/stds(ii));
            end
        case {'logmean'}
            if isfield(params,'nonLinOutMean')
                means = params.nonLinOutMean;
            else
                means = nanmean(S);
                params.nonLinOutMean = means;
            end
            % For logmean, d is also a small value to avoid -Inf
            d = params.nonLinOutParam;
            for ii=1:size(S,2)
                S(:,ii) = log(d+S(:,ii)/means(ii));
            end
        case {'linear'}
            disp('Linearity requested. Do nothing.');
    end
else
    if length(params.nonLinOutExp) == 1
        S = abs(S).^params.nonLinOutExp.*sign(S);
    else
        Pt = [];
        for ii=1:length(params.nonLinOutExp)
            Pt = cat(2, Pt, abs(S).^params.nonLinOutExp(ii).*sign(S));
        end
        S = Pt;
    end
end

if isfield(params,'verbose') && params.verbose
    fprintf(' done.\n'); 
end

params.nChan = size(S,2);

% Output
varargout{1} = S;
if nargout>1
    varargout{2} = params;
end
