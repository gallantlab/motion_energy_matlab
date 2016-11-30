function varargout = preprocWavelets_grid(S, params)
% Usage: [Spreproc, params] = preprocWavelets_grid(S, params);
%
% A script for preprocessing of stimuli using a Gabor wavelet basis set
%
% INPUT:
%                 [S] = A X-by-Y-by-T matrix containing stimuli (movie);
%                       range of images will affect filter output values
%                       (0-1 vs 0-255), but that difference will go away
%                       after normalization of channels (if performed in
%                       subsequent preprocessing steps)
%            [params] = structure that contains parameters for
%                           preprocessing, with fields: 
%  --- STRFlab housekeeping ---
%    .show_or_preprocess = If this is set to 0, the function returns wavelets
%                       of size(S) * number of channels, instead of preprocessed
%                       stimuli. This may be used for visualization purpose.
%                       If .valid_w_index is also set, this returns only a subset of
%                       wavelets specified by .valid_w_index. (default: 1)
%  --- Misc ---
% 
%         .f_step_log = A flag to specify linear (false) or log (true) step of frequency (default: 0)
%      .valid_w_index = This is used to specify a subset of wavelets to obtain.
%                       (See .show_or_preprocess)
%            .verbose = a flag for verbose mode (default: 1)
%           .zeromean = a flag for whether to make stimulus zero mean
%                       before preprocessing (default: 0)
%           fenv_mode = a flag to set fenv_mode (no diff btw. spatial and
%                       temporal envelopes??) (default: 0)
%      gaborcachemode = (??) (default: 0)
%
%  --- Orientation / direction parameters ---
%       .dirdivisions = Number of directions for wavelets (default: 8)
% .directionSelective = a flag for whether the model should be selective
%                       for different directions of motion at the same orientation
% 
%  --- Spatial freqeuncy / location parameters ---
%        .sfdivisions = Number of spatial frequencies (default: 5)
%              .sfmax = The maximum spatial frequency/stimulus size at zero velocity (default: 9)
%              .sfmin = The minimum spatial frequency/stimulus size  at zero velocity (default: 2)
%           .local_dc = A flag to add localized dc (i.e., 0 spatial freq.) channels (default: 0)
%      .sf_gaussratio = The ratio between the Gaussian window and spatial frequency (default: 0.5)
%           .std_step = Spatial separation of each wavelet in terms of sigma of
%                       the Gaussian window (default: 2.5)
%           .senv_max = The maximum spatial envelope (default: 0.3)
%           .wrap_all = Whether to cluster channels tightly at the center
%                       of the screen (1) or not (0). (default: 0)
%  --- Temporal freqeuncy parameters ---
%        .tfdivisions = Number of velocities (default: 5)
%              .tsize = Number of frames to calculate wavelets (default: 10)
%              .tfmax = The maximum temporal frequency/stimulus size (default: 3.5)
%      .tf_gaussratio = The ratio between the Gaussian window and temporal frequency (default: 0.4)
%           .tenv_max = The maximum temporal envelope (default: 0.3)
%
%  --- Phase parameters ---
%          .phasemode = A parameter to specify how to deal with phase information 
%                       (i.e., how to combine [or not] Gabor channels with different phases)
%                         0: spectral amplitude (default)
%                         1: linear sin and cos phase ampliture (2x number of wavelets)
%                         2: half-rectified sin and cos phase amplitude (4x number of wavelets)
%                         3: 0+1 (3x number of wavelets)
%                         4: 0+2 (5x number of wavelets)
%                         5: phase, atan2(chout90,chout0)
%                         6: 0+5 (2x number of wavelets)
%                         7: phase, atan2(chout90,chout0), half-rectified
%                         8: 0+7 (3x number of wavelets)
%    .phasemode_sfmax = The maximum spatial frequency to use phase information
%                       For higher frequency wavelets over this value, only
%                       spectral amplitude information are used.  (default: Inf)
%
% OUTPUT:
%          [Spreproc] = Preprocessed stimuli that can be used for STRF fitting.
%                       NxD matrix, N=sample size, D=dimensions (model channels)
%            [params] = structure that contains parameters for preprocessing, with additional fields:
%              .nChan = Number of preprocessed channels
%                       (dimensionality of each data vector, AKA: D=dimensions)
%        .gaborparams = A set of parameters for each Gabor wavelet.
%                       This is a p-by-D matrix where p is number of parameters (8)
%                       and D is the number of wavelet channels
%                       Each field in gaborparams represents:
%                       [pos_x pos_y direction s_freq t_freq s_size t_size phasevalue]
%                       phasevalue can be 0 to 6, where
%                         0: spectra
%                         1: linear sin transform
%                         2: linear cos transform
%                         3: half-rectified sin transform (positive values)
%                         4: half-rectified sin transform (negative values)
%                         5: half-rectified cos transform (positive values)
%                         6: half-rectified cos transform (negative values)
%                         7: dPhase/dt
%                         8: dPhase/dt (positive values)
%                         9: dPhase/dt (negative values)
%     .zeromean_value = Offset value of total movie. This is set if .zeromean is non-zero.
%
% EXAMPLE:
%  params = preprocWavelets;
%    returns the default set of parameters.
%
%  [Spreproc, params] = preprocWavelets_grid(S, PARAMS)
%    returns preprocessed stimuli (or wavelets) and parameters
%
% SEE ALSO: make3dgabor, preprocSpectra
% ====================


% Default parameters
% STRFlab housekeeping
dParams.show_or_preprocess = 1;
% Misc
dParams.f_step_log = 0;
dParams.fenv_mode = 0; %??
dParams.gaborcachemode = 0;
dParams.valid_w_index = NaN;
dParams.zeromean = 1;
dParams.verbose = 1;
% Orientation / direction 
dParams.dirdivisions = 8;
dParams.directionSelective = 1;
% Spatial frequency / location
dParams.sfdivisions = 5;
dParams.sfmax = 9.0;
dParams.sfmin = 2.0;
dParams.local_dc = 0;
dParams.sf_gaussratio = 0.5;
dParams.senv_max = 0.3;
dParams.std_step = 2.5;
dParams.wrap_all = 0;
% Temporal frequency
dParams.tfmax = 3.0;
dParams.tfmin = 1.0;
dParams.tfdivisions = 5;
dParams.tsize = 9;
dParams.zerotf = 1;
dParams.tf_gaussratio = 0.4;
dParams.tenv_max = 0.3;
% Phase
dParams.phasemode = 0;
dParams.phasemode_sfmax = Inf; % linked to .phasemode
% Fill in default parameters
if ~exist('params','var')
    params = struct;
end
params = defaultOpt(params,dParams);
params.class = 'preprocWavelets_grid';
% Misc fill-in
if isfield(params,'fenv_mode')
    if ~isfield(params, 'f_gaussratio')
        params.f_gaussratio = 0.5;
    end
    if ~isfield(params, 'fenv_max')
        params.fenv_max = 0.3;
    end
end
% Return parameters if no input
if nargin<1 
    varargout{1} = params;
    return;
end

% Timing
start_t = cputime;
% Stimulus check
stimxytsize = size(S);
% Stimulus aspect ratio; always X/Y
aspect_ratio = stimxytsize(2)/stimxytsize(1);
if length(stimxytsize) == 2
    stimxytsize = [stimxytsize 1]; % make sure 3 dimensions
elseif length(stimxytsize)==4
    if params.verbose
        disp('Processing color channels separately...')
    end
    nColorChannels = size(S,4);
    % Color stimuli; recursive call to process color channels separately
    Spreproc = [];
    for ci=1:nColorChannels
        [tstim p] = feval(params.class, S(:,:,:,ci), params);
        Spreproc = cat(2,Spreproc, tstim);
    end
    params = p;
    varargout{1} = Spreproc;
    if nargout >1
        varargout{2} = params;
    end
    % Done!
    return
end

patchxytsize = [stimxytsize(1:2) params.tsize];
xypixels = prod(patchxytsize(1:2));
verbose = params.verbose;

if ~isfloat(S)
    S = single(S);
end
S = reshape(S, [prod(stimxytsize(1:2)) stimxytsize(3)]);

if params.show_or_preprocess
    if params.zeromean
        if verbose, fprintf('[[zero mean stimuli]]\n'); end
        if isfield(params, 'zeromean_value')
            S = bsxfun(@minus, S, params.zeromean_value);
        else
            thismean = mean(S(:));
            S = bsxfun(@minus, S, thismean);
            params.zeromean_value = thismean;
        end
    end
end

% Make a list of gabor parameters
if ~isfield(params,'gaborparams') || params.phasemode == 5 || ...
        params.phasemode == 7
    if verbose, fprintf('Making a list of gabor parameters... '); end
    % Added aspect ratio as necessary influence on Gabor parameters
    [gparams] = getGaborParameters(params,aspect_ratio);
else
    gparams = params.gaborparams;
end
waveletchannelnum = size(gparams,2);

if verbose, fprintf('channel num: %d\n', waveletchannelnum); end

if verbose && any(params.valid_w_index)
    fprintf('Valid channel num: %d\n', length(params.valid_w_index));
end


% Set up a matrix to fill in
if params.show_or_preprocess
    if verbose, disp('Preprocessing...'); end
    Spreproc = zeros(stimxytsize(3), waveletchannelnum, 'single');
else
    if verbose, disp('Making wavelets...'); end
    if ~any(params.valid_w_index)
        gnum = length(waveletchannelnum);
    else
        gnum = length(params.valid_w_index);
    end
    gaborbank = zeros([patchxytsize gnum], 'single');
end

%---------------------------------------------------------------------
% Preprocessing
%---------------------------------------------------------------------
% ignore wavelet pixels for speed-up where:
masklimit = 0.001;   %% pixel value < masklimit AND
maskenv_below = 0.1; % spatial envelope < maskenv_below x stimulus size

if params.gaborcachemode==1
    gaborcache = zeros([2 prod(patchxytsize(1:2)) waveletchannelnum], 'single');
    gtwcache = zeros([2 params.tsize waveletchannelnum], 'single');
end

lastgparam = zeros(9,1);
wcount = 0;
for ii=1:waveletchannelnum
    
    if any(params.valid_w_index) && ~any(ii==params.valid_w_index), continue, end
    
    thisgparam = gparams(:,ii);
    thesame = 1;
    if any(thisgparam([1:7,9]) ~=lastgparam([1:7,9])) 
        thesame = 0;
    end
    if ~thesame
        if params.gaborcachemode==2
            gabors = params.gaborcache(:,:,ii);
            gtw = params.gtwcache(:,:,ii);
        else
            [gabor0,gabor90,gtw] = make3dgabor_frames(patchxytsize, [thisgparam(1:7); 0; thisgparam(9)]);
            gabors = [gabor0(:) gabor90(:)]';
        end
        if params.gaborcachemode==1
            gaborcache(:,:,ii) = gabors;
            gtwcache(:,:,ii) = gtw;
        end
        lastgparam = thisgparam;
    end
    phaseparam = thisgparam(8);
    if params.show_or_preprocess
        if ~thesame
            senv = thisgparam(6);
            if senv<maskenv_below
                smask = find(sum(abs(gabors),1)>masklimit);
                [chout0,chout90] = dotdelay_frames(gabors(:,smask), gtw, S(smask,:));
            else
                [chout0,chout90] = dotdelay_frames(gabors, gtw, S);
            end
        end
        switch phaseparam
            case 0
                chout = sqrt(chout0.^2 + chout90.^2);
                Spreproc(:,ii) = chout;
            case 1
                chout = chout0;
                Spreproc(:,ii) = chout;
            case 2
                chout = chout90;
                Spreproc(:,ii) = chout;
            case 3
                chout = chout0;
                chout(chout<0) = 0;
                Spreproc(:,ii) = chout;
            case 4
                chout = chout0;
                chout(chout>0) = 0;
                Spreproc(:,ii) = -chout;
            case 5
                chout = chout90;
                chout(chout<0) = 0;
                Spreproc(:,ii) = chout;
            case 6
                chout = chout90;
                chout(chout>0) = 0;
                Spreproc(:,ii) = -chout;
            case 7
                chout = atan2(chout90,chout0);
                dtphase = [0; diff(chout,1,1)];
                dtphase = dtphase+ -2*pi*sign(dtphase).*round(abs(dtphase)./(2*pi));
                Spreproc(:,ii) = dtphase;
            case 8
                chout = atan2(chout90,chout0);
                dtphase = [0; diff(chout,1,1)];
                dtphase = dtphase+ -2*pi*sign(dtphase).* ...
                    round(abs(dtphase)./(2*pi));
                dtphase(dtphase<0) = 0;
                Spreproc(:,ii) = dtphase;
            case 9
                chout = atan2(chout90,chout0);
                dtphase = [0; diff(chout,1,1)];
                dtphase = dtphase+ -2*pi*sign(dtphase).* ...
                    round(abs(dtphase)./(2*pi));
                dtphase(dtphase>0) = 0;
                Spreproc(:,ii) = -dtphase;
        end
    else
        wcount = wcount + 1;
        switch phaseparam
            case {1,3,4}
                % reconstruct space-time Gabor
                rgs = gabors(1,:)'*gtw(2,:)+gabors(2,:)'*gtw(1,:);
                rgs=reshape(rgs, [patchxytsize]);
                gaborbank(:,:,:,wcount) = rgs;
            case {0,2,5,6,7,8,9}
                % reconstruct space-time Gabor
                rgc = -gabors(1,:)'*gtw(1,:)+gabors(2,:)'*gtw(2,:);
                rgc=reshape(rgc, [patchxytsize]);
                gaborbank(:,:,:,wcount) = rgc;
        end
    end
    
    if verbose
        progressdot(ii,50,1000,waveletchannelnum);
    end
end

if params.gaborcachemode==1
    params.gaborcache = gaborcache;
    params.gtwcache = gtwcache;
    params.gaborcachemode = 2;
end

if verbose
    disp(sprintf('Wavelet preprocessing done in %.1f min (cputime).', (cputime-start_t)/60));
    if params.show_or_preprocess
        disp(sprintf('%d channels, %d samples', size(Spreproc,2), size(Spreproc,1)));
    else
        disp(sprintf('%d channels', size(gaborbank,4)));
    end
end

if params.show_or_preprocess
    if params.phasemode==5 || params.phasemode==6 || params.phasemode==7 ...
            || params.phasemode==8
        pind = find(gparams(8,:)==7 | gparams(8,:)==8 | gparams(8,:)==9);
        disp('thresholding phase channels...');
        for p=1:length(pind)
            phasech = Spreproc(:,pind(p));
            if gparams(8,pind(p)-1) == 0 % look for the
                % corresponding amplitude channel
                ampch = Spreproc(:,pind(p)-1);
            else
                ampch = Spreproc(:,pind(p)-2);
            end
            a_thresh = nanstd(ampch)*params.a_thresh;
            avalind = ampch>a_thresh;
            avalind = and(avalind, [0; avalind(1:end-1)]);
            phasech(~avalind) = 0;
            Spreproc(:,pind(p)) = phasech;
        end
        if params.phasemode==5 || params.phasemode==7 % return dPhase/dt channels only
            Spreproc = Spreproc(:,pind);
            gparams = gparams(:,pind);
            fprintf('Using only dPhase/dt channels: %d\n', size(Spreproc,2));
        end
    end
    
else % return gabors, not pre-processed data
    Spreproc = gaborbank;
end

params.gaborparams = gparams;
params.nChan = size(gparams,2);

varargout{1} = Spreproc;
if nargout>1
    varargout{2} = params;
end

return;


%---------------------------------------------------------------------
% Making a list of gabor parameters
%---------------------------------------------------------------------
function gparams = getGaborParameters(params,aspect_ratio)


if params.f_step_log
    sf_array = logspace(log10(params.sfmin), log10(params.sfmax), params.sfdivisions);
    if params.zerotf
        tf_array = logspace(log10(params.tfmin), log10(params.tfmax), params.tfdivisions-1);
        tf_array = [0 tf_array];
    else
        tf_array = logspace(log10(params.tfmin), log10(params.tfmax), params.tfdivisions);
    end
else
    sf_array = linspace(params.sfmin, params.sfmax, params.sfdivisions);
    tf_array = linspace(params.tfmin, params.tfmax, params.tfdivisions);
    
end

dir_array = (0:params.dirdivisions-1)/params.dirdivisions*360;

switch params.phasemode
    case 0  % spectral amplitudes
        pmarray = [0];
    case 1  % linear sin and cos transform amplitudes
        pmarray = [1 2];
    case 2  % half rectified sin and cos amplitudes
        pmarray = [3 4 5 6];
    case 3  % 0+1
        pmarray = [0 1 2];
    case 4  % 0+2
        pmarray = [0 3 4 5 6];
    case 5  % phase: atan2(sin,cos)
        pmarray = [0 7];
    case 6  % 0+5
        pmarray = [0 7];
    case 7  % phase: atan2(sin,cos), half-rectified
        pmarray = [0 8 9];
    case 8  % 0+7
        pmarray = [0 8 9];
end

dirstart = 1;
if params.local_dc
    dirstart = 0; %% add local dc channels
end

waveletcount = 0;
gparams = zeros(8, 20000, 'single'); % prepare for some amount of memory for gparams
% Add a row to gparams to account for aspect ratio
gparams = [gparams;ones(1,size(gparams,2),'single')];

for ti=1:params.tfdivisions
    tf = tf_array(ti);
    for fi = 1:params.sfdivisions
        sf = sf_array(fi);
        
        if params.fenv_mode
            f = sqrt(sf.^2+tf.^2);
            fenv = min([params.fenv_max 1/f*params.f_gaussratio]);
            tenv = fenv;
            senv = fenv;
        else
            senv = params.senv_max;
            if sf ~= 0
                senv = min([params.senv_max 1/sf*params.sf_gaussratio]);
            end
            tenv = params.tenv_max;
            if tf ~= 0
                tenv = min([params.tenv_max 1/tf*params.tf_gaussratio]);
            end
        end
        
        if params.directionSelective == 0
            tf = tf + i;
        end
        
        % Account for asymmetrical images
        if aspect_ratio==1
            % Symmetrical images
            numsps2 =floor((1-senv*params.std_step)/(params.std_step*senv)/2);
            numsps2 = max([numsps2 0]);
            if numsps2>=1 && params.wrap_all
                numsps2 = numsps2 + 1;
            end
            centers = senv*params.std_step*(-numsps2:numsps2) + 0.5;
            [cx, cy] = meshgrid(centers, centers);    
        else
            % aspect_ratio is x/y. Thus ar*x = true x OR y/ar = true y
            % Compute 
            g_sz_x = senv*params.std_step;
            n_gabors_x = floor((1-g_sz_x)/(g_sz_x)/2);
            n_gabors_x = max([n_gabors_x,0]);
            g_sz_y = senv*params.std_step*aspect_ratio;
            n_gabors_y = floor((1-g_sz_y)/(g_sz_y)/2);
            n_gabors_y = max([n_gabors_y,0]);
            %REPLACED:
            %numsps2 =floor((1-senv*params.std_step)/(params.std_step*senv)/2);
            %numsps2 = max([numsps2 0]);
            %if numsps2>=1 && params.wrap_all
            if params.wrap_all
                error('I don''t know what to do with wrap_all parameter yet w/ asymmetrical images...')
                %numsps2 = numsps2 + 1;
            end
            centers_x = g_sz_x*(-n_gabors_x:n_gabors_x) + 0.5;
            centers_y = g_sz_y*(-n_gabors_y:n_gabors_y) + 0.5;
            [cx, cy] = meshgrid(centers_x, centers_y);
            % Elongate gabors if differential sampling in x and y does not
            % make the gabors circular
            % BUT ONLY IF there is no elongation parameter? (BUT WHERE THE
            % HELL DOES/DID THAT COME IN? nowhere that I (ML) can find in
            % Shinji's code. Must have been a hand-coded addition to
            % make3dgabor_frames.
            sampling_aspect_ratio = length(centers_x)/length(centers_y);
            elong = aspect_ratio / sampling_aspect_ratio;
            %keyboard;
        end        
        thisnumdirs = length(dir_array);
        if tf == 0 || params.directionSelective == 0
            thisnumdirs = ceil(thisnumdirs/2);  % use only ~180 deg
        end
        if sf == 0
            thisnumdirs = 1;
        end
        for xyi = 1:length(cx(:))
            xcenter = cx(xyi);
            ycenter = cy(xyi);
            for diri = dirstart:thisnumdirs
                if diri
                    dir = dir_array(diri); thissf = sf;
                else
                    if params.local_dc == 1
                        dir = 0; thissf = 0; % local dc channels
                    else
                        dir = 0; thissf = sf*0.01; % to avoid the exact same channel
                    end
                end
                if  thissf >= params.phasemode_sfmax
                    waveletcount = waveletcount+1;
                    thisgparam = [xcenter ycenter dir thissf tf senv tenv 0 1];
                    if aspect_ratio~=1
                        thisgparam(9) = max(elong,1);
                    end
                    gparams(:,waveletcount) = thisgparam;
                else
                    for pmod = pmarray
                        waveletcount = waveletcount+1;
                        thisgparam = [xcenter ycenter dir thissf tf senv tenv pmod 1];
                        if aspect_ratio~=1
                            thisgparam(9) = max(elong,1);
                        end
                        gparams(:,waveletcount) = thisgparam;
                    end
                end
            end
        end
    end
end

gparams = gparams(:,1:waveletcount);

