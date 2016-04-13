function params = preprocColorSpace_GetMetaParams(argNum)

params.class = 'preprocColorSpace';

switch argNum
    case 1
        % Default grayscale space
        params.class = 'preprocColorSpace';
        params.colorconv = 'LAB<-RGB';
        params.colorchannels = 1;
        params.gamma = 1.0;
        params.verbose = true;
    otherwise
        error('Unknown argument!')
end