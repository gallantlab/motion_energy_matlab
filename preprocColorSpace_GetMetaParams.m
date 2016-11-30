function params = preprocColorSpace_GetMetaParams(argNum)

params.class = 'preprocColorSpace';

switch argNum
    case 1
        % Convert to L*A*B colorspace, keep luminance channel
        params.class = 'preprocColorSpace';
        params.colorconv = 'rgb2lab';
        params.colorchannels = 1;
        params.gamma = 1.0;
        params.verbose = true;
    otherwise
        error('Unknown argument!')
end