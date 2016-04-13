function LMS_out = rgb_to_lms(rgb_im, t, g)

    origsz = size(rgb_im);
    if ndims(rgb_im) == 3
        rgb_im = reshape(rgb_im, [], 3);
    end
    rgb_im = double(rgb_im) + 1;
    gr = g(:,1)';
    gg = g(:,2)';
    gb = g(:,3)';

    % modified by SN: assuming gamma corrected in advance
    %rgb_gamma = [gr([rgb_im(:,1)])' gg([rgb_im(:,2)])' gb([rgb_im(:,3)])'];
    rgb_gamma = rgb_im;
    
    % rgb_gamma = [gamma([rgb_im(:,1)]) gamma([rgb_im(:,2)])
    % gamma([rgb_im(:,3)])];
    % E = RGB(:,2)*rgb_gamma(:,1)' + RGB(:,3)*rgb_gamma(:,2)' +
    % RGB(:,4)*rgb_gamma(:,3)';
    % LMS_out = E' * LMS(:,2:4);
    LMS_out = rgb_gamma * t;

    LMS_out = reshape(LMS_out, origsz);