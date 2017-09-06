function params = preprocWavelets_grid_GetMetaParams(Arg)
% Usage: params = preprocWavelets_grid_GetMetaParams(Arg)
% 
% Returns parameter set for a preprocWavelets_grid. The parameters
% specified by setting Arg to 2 are the parameters used to compute
% motion energy in Nishimoto et al (2011); Arg = 1 is a similar model
% with fewer channels that provides comparable results for modeling
% fMRI data.
%
% See code for other parameter presets. 

params = preprocWavelets_grid;
params.argNum = Arg;

switch Arg
    case 1
        % Motion energy model with fewer features
        % Housekeeping
        params.class = 'preprocWavelets_grid';
        params.show_or_preprocess = 1; % True to preprocess; false to return gabor channels
        params.verbose = 1;
        params.gaborcachemode = 0; % whether or not to cache calculated results (could become faster)
        params.valid_w_index = NaN; % Select particular gabor channels by number
        % Temporal frequency params
        params.tfdivisions = 3; % Number of temporal frequencies; [tfmin...tfmax] or [0, tfmin...tfmax] if zerorf=1
        params.tfmax = 2.66667; % = 4hz @ 15 fps ([tfmax] cycles per [tsize] frames at 15 fps; 2.66667/10*15 = 4 Hz) 
        params.tfmin = 1.33333; % = 2hz @ 15 fps (1.33333/10*15 = 2 Hz)
        params.tsize = 10; % The size of temporal window (frames)
        params.tf_gaussratio = 10; % temporal frequency to gauss envelope ratio of Gabor; bigger number = more waves (larger envelope)
        params.tenv_max = 0.3000; % the maximum gaussian envelope size (relative to tsize)
        params.zerotf = 1; % Include 0 Hz (static) energy channels
        params.f_gaussratio = .5; % frequency to gauss ratio of Gabor; obsolete
        % Orientation/direction params
        params.dirdivisions = 8;
        params.local_dc = 1; % T/F. Include circular gaussians (w/ no spat. freq.)
        params.directionSelective = 1;
        % Spatial extent params
        params.sfdivisions = 5;
        params.sfmax = 24; %
        params.sfmin = 1.5; %
        params.f_step_log = 1; % Applies to both SF and TF?
        params.std_step = 4; % Governs how closely spaced channels are; a reasonable range is 2.5-4
        params.sf_gaussratio = 0.6000; % 81 channels @maxsf=24; 9x9 ; 13x13 @maxsf=32
        params.fenv_mode = 0; % use same env for spatial & temporal gabors
        params.senv_max = 0.3000;
        params.wrap_all = 0; % whether or not the filters cover the very edge of images
        % Handling phase
        params.phasemode = 0; % Determines how to do phase (square & sum quadrature pairs (i.e., energy), linear, rectified, etc)
        params.phasemode_sfmax = NaN; % Calculate only energy channels (e.g., no linear channels) if sf exceeds this number.
        params.zeromean = 1;
    case 2
        % larger motion energy model 
        % (as in Nishimoto et al 2011)
        % STRFlab conventions, housekeeping
        params.class = 'preprocWavelets_grid';
        params.show_or_preprocess = 1;
        params.wrap_all = 0;
        params.verbose = 1;
        params.gaborcachemode = 0;
        params.valid_w_index = NaN;
        % Temporal frequency params
        params.tfdivisions = 3;
        params.tfmax = 2.66667;
        params.tfmin = 1.33333;
        params.tsize = 10;
        params.tf_gaussratio = 10; 
        params.tenv_max = 0.3000;
        params.zerotf = 1;
        % Orientation/direction params
        params.dirdivisions = 8;
        params.local_dc = 1; 
        params.directionSelective = 1;
        % Spatial extent params
        params.sfdivisions = 5;
        params.sfmax = 32; %
        params.sfmin = 2; %
        params.f_step_log = 1; % Applies to both SF and TF?
        params.std_step = 3.5; % Governs how closely spaced channels are
        params.sf_gaussratio = 0.6000; % 81 channels @maxsf=24; 9x9 ; 13x13 @maxsf=32
        params.fenv_mode = 0; % (whether to use fenv_max for both senv_max and tenv_max) 
        params.senv_max = 0.3000;
        % Nonlinearities
        params.phasemode = 0; % Determines how to do phase (square & sum quadrature pairs, etc)
        params.phasemode_sfmax = NaN; % No idea
        params.zeromean = 1;        
        
    case 3
        % Same as 1, but NO PYRAMID (high spatial frequency channels only)
        % smaller motion energy model
        % STRFlab conventions, housekeeping
        params.wclass = 'preprocWavelets_grid';
        params.show_or_preprocess = 1;
        params.wrap_all = 0;
        params.verbose = 1;
        params.gaborcachemode = 0;
        % Temporal frequency params
        params.tfdivisions = 3;
        params.tfmax = 2.66667;
        params.tfmin = 1.33333;
        params.tsize = 10;
        params.tf_gaussratio = 10; 
        params.tenv_max = 0.3000;
        params.zerotf = 1;
        params.f_gaussratio = .5;
        % Orientation/direction params
        params.dirdivisions = 8;
        params.local_dc = 1; % T/F. Include gaussians w/ no spat. freq.
        params.directionSelective = 1;
        % Spatial extent params
        params.sfdivisions = 1;
        params.sfmax = 24; %
        params.sfmin = 24; %
        params.f_step_log = 1; % Applies to both SF and TF?
        params.std_step = 4; % Governs how closely spaced channels are
        params.sf_gaussratio = 0.6000; % 81 channels @maxsf=24; 9x9 ; 13x13 @maxsf=32
        params.fenv_mode = 0; % no idea
        params.senv_max = 0.3000;
        % Nonlinearities
        params.phasemode = 0; % Determines how to do phase (square & sum quadrature pairs, etc)
        params.phasemode_sfmax = NaN; % No idea
        params.zeromean = 1;
    case 4
        % Same as 2, but NO PYRAMID (high spatial frequency channels only)
        % STRFlab conventions, housekeeping
        params.class = 'preprocWavelets_grid';
        params.show_or_preprocess = 1;
        params.wrap_all = 0;
        params.verbose = 1;
        params.gaborcachemode = 0;
        params.valid_w_index = NaN;
        % Temporal frequency params
        params.tfdivisions = 3;
        params.tfmax = 2.66667;
        params.tfmin = 1.33333;
        params.tsize = 10;
        params.tf_gaussratio = 10; 
        params.tenv_max = 0.3000;
        params.zerotf = 1;
        % Orientation/direction params
        params.dirdivisions = 8;
        params.local_dc = 1; 
        params.directionSelective = 1;
        % Spatial extent params
        params.sfdivisions = 1;
        params.sfmax = 32; %
        params.sfmin = 32; %
        params.f_step_log = 1; % Applies to both SF and TF?
        params.std_step = 3.5; % Governs how closely spaced channels are
        params.sf_gaussratio = 0.6000; % 81 channels @maxsf=24; 9x9 ; 13x13 @maxsf=32
        params.fenv_mode = 0; % (whether to use fenv_max for both senv_max and tenv_max) 
        params.senv_max = 0.3000;
        % Nonlinearities
        params.phasemode = 0; % Determines how to do phase (square & sum quadrature pairs, etc)
        params.phasemode_sfmax = NaN; % max (min?) SF for which phase is computed
        params.zeromean = 1; 
    case 5
        % same as 1 w/ NO TEMPORAL CHANNELS (static Gabor wavelet model)
        % STRFlab conventions, housekeeping
        params.class = 'preprocWavelets_grid';
        params.show_or_preprocess = 1; % True to preprocess; false to return gabor channels
        params.verbose = 1;
        params.gaborcachemode = 0;
        params.valid_w_index = NaN; % Select particular gabor channels by number
        % Temporal frequency params
        params.tfdivisions = 1;
        params.tfmax = 0; %2.66667; % = 4hz @ 15 fps
        params.tfmin = 0; %1.33333; % = 2hz @ 15 fps
        params.tsize = 1;
        params.tf_gaussratio = 1; %10; 
        params.tenv_max = 0.3000;
        params.zerotf = 1;
        params.f_gaussratio = .5;
        % Orientation/direction params
        params.dirdivisions = 8;
        params.local_dc = 1; % T/F. Include circular gaussians (w/ no spat. freq.)
        params.directionSelective = 1;
        % Spatial extent params
        params.sfdivisions = 5;
        params.sfmax = 24; %
        params.sfmin = 1.5; %
        params.f_step_log = 1; % Applies to both SF and TF?
        params.std_step = 4; % Governs how closely spaced channels are; a reasonable range is 2.5-4
        params.sf_gaussratio = 0.6000; % 81 channels @maxsf=24; 9x9 ; 13x13 @maxsf=32
        params.fenv_mode = 0; % use same env for spatial & temporal gabors
        params.senv_max = 0.3000;
        params.wrap_all = 0;
        % Handling phase
        params.phasemode = 0; % Determines how to do phase (square & sum quadrature pairs, etc)
        params.phasemode_sfmax = NaN; % No idea
        params.zeromean = 1;
    case 6
        % same as 2 w/ NO TEMPORAL CHANNELS (static Gabor wavelet model)
        % (Approximately as in Kay et al 2008, Naselaris et al 2009)
        % STRFlab conventions, housekeeping
        params.class = 'preprocWavelets_grid';
        params.show_or_preprocess = 1;
        params.wrap_all = 0;
        params.verbose = 1;
        params.gaborcachemode = 0;
        params.valid_w_index = NaN;
        % Temporal frequency params
        params.tfdivisions = 1;
        params.tfmax = 0;
        params.tfmin = 0;
        params.tsize = 1;
        params.tf_gaussratio = 1; 
        params.tenv_max = 0.3000;
        params.zerotf = 1;
        % Orientation/direction params
        params.dirdivisions = 8;
        params.local_dc = 1; 
        params.directionSelective = 1;
        % Spatial extent params
        params.sfdivisions = 5;
        params.sfmax = 32; %
        params.sfmin = 2; %
        params.f_step_log = 1; % Applies to both SF and TF?
        params.std_step = 3.5; % Governs how closely spaced channels are
        params.sf_gaussratio = 0.6000; % 81 channels @maxsf=24; 9x9 ; 13x13 @maxsf=32
        params.fenv_mode = 0; % (whether to use fenv_max for both senv_max and tenv_max) 
        params.senv_max = 0.3000;
        % Nonlinearities
        params.phasemode = 0; % Determines how to do phase (square & sum quadrature pairs, etc)
        params.phasemode_sfmax = NaN; % No idea
        params.zeromean = 1;        
    case 7
        % Motion energy model for HCP data
        % STRFlab conventions, housekeeping
        pp.class = 'preprocWavelets_grid';
        pp.show_or_preprocess = 1;
        pp.wrap_all = 0;
        pp.verbose = 1;
        pp.gaborcachemode = 0;
        pp.valid_w_index = NaN;
        % Temporal frequency params
        pp.tfdivisions = 3;
        pp.tfmax = 1.66667; % = 4hz @ 24 fps, w/ 10 frame t limit
        pp.tfmin = 0.83333; % = 2hz @ 24 fps, w/ 10 frame t limit
        pp.tsize = 10;
        pp.tf_gaussratio = 10; 
        pp.tenv_max = 0.3000;
        pp.zerotf = 1;
        % Orientation/direction params
        pp.dirdivisions = 8;
        pp.local_dc = 1; 
        pp.directionSelective = 1;
        % Spatial extent params
        pp.sfdivisions = 5;
        pp.sfmax = 32; %
        pp.sfmin = 2; %
        pp.f_step_log = 1; % Applies to both SF and TF?
        pp.std_step = 3.5; % Governs how closely spaced channels are
        pp.sf_gaussratio = 0.6000; % 81 channels @maxsf=24; 9x9 ; 13x13 @maxsf=32
        pp.fenv_mode = 0; % (whether to use fenv_max for both senv_max and tenv_max) 
        pp.senv_max = 0.3000;
        % Nonlinearities
        pp.phasemode = 0; % Determines how to do phase (square & sum quadrature pairs, etc)
        pp.phasemode_sfmax = NaN; % No idea
        pp.zeromean = 1;        
    otherwise
        error('Unknown argument!')        
end
