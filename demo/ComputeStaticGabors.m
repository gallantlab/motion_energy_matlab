% Load images

d = load('NaselarisStansburyGallant_2012_ValImages.mat');
% the field d.S is an array that is (128 x 128 x 3 x 126); (X x Y x Color x
% Images).  The images are stored as 8-bit integer arrays (no decimal
% places, with pixel values from 0-255). These should be converted to
% floating point decimals from 0-1:
S  = single(d.S)/255;

%% Preprocessing
% Conver to grayscale (luminance only)
% The argument 1 here indicates a pre-specified set of parameters to feed
% to thie preprocColorSpace function to convert from RGB images to 
% luminance values by converting from RGB to L*A*B colorspace and then
% keeping only the luminance channel. (You could also use matlab's
% rgb2gray.m function, but this is more principled.) Inspect cparams to see
% what those parameters are.
cparams = preprocColorSpace_GetMetaParams(1);
[S_lum, cparams] = preprocColorSpace(S, cparams);

%% Gabor wavelet processing
% Process with Gabor wavelets
% The numerical argument here specifies a set of parameters for the
% preprocWavelets_grid function, that dictate the locations, spatial
% frequencies, phases, and orientations of Gabors to use. 6 specifies Gabor
% wavelets with no temporal component (suitable for static images). 
gparams = preprocWavelets_grid_GetMetaParams(6);
[S_gab, gparams] = preprocWavelets_grid(S_lum, gparams);

%% Optional additions
% Compute log of each channel to scale down very large values
nlparams = preprocNonLinearOut_GetMetaParams(1);
[S_nl, nlparams] = preprocNonLinearOut(S_gab, nlparams);

% Z-score each channel
nrmparams = preprocNormalize_GetMetaParams(1);
[S_fin, nrmparams] = preprocNormalize(S_nl, nrmparams);

%% Display output
disp('Final matrix size (images x features):')
disp(size(S_fin));

imagesc(S_fin);
ylabel('Image')
xlabel('Gabor wavelet feature')
caxis([-3,3]);
colorbar();