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
    case 2
        % Convert to grayscale using rgb2gray (inferior, but present on
        % older matlab versions
        params.class = 'preprocColorSpace';
        params.colorconv = 'rgb2gray';
        params.colorchannels = 1;
        params.gamma = 1.0;
        params.verbose = true;        
    otherwise
        error('Unknown argument!')
end