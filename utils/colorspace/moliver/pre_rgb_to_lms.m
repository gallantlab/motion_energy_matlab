function [t, gamma] = pre_rgb_to_lms(model)
    
    if model == 1
        load('/auto/k1/moliver/code/color_strf/RGB.mat')
        load('/auto/k1/moliver/code/color_strf/LMS.mat')
        load('/auto/k1/moliver/code/color_strf/avo.mat')
        t = RGB(:,2:4)'*LMS(:,2:4);

        xx = [x; x];
        rr = [r1; r2];
        gg = [g1; g2];
        bb = [b1; b2];
        rgb = [rgb1; rgb2];

        xall = [0:255];

        [pr, sr] = polyfit(xx,rr,2);
        [pg, sg] = polyfit(xx,gg,2);
        [pb, sb] = polyfit(xx,bb,3);

        yr = polyval(pr, xall);
        yr = yr./max(yr);
        yg = polyval(pg, xall);
        yg = yg./max(yg);
        yb = polyval(pb, xall);
        yb = yb./max(yb);
        gamma = [yr' yg' yb'];
        gamma(find(gamma < 0)) = 0.0001;
    elseif model == 2
        load('/auto/k1/moliver/code/color_strf/RGB.mat')
        load('/auto/k1/moliver/code/color_strf/LMS.mat')
        t = RGB(:,2:4)'*LMS(:,2:4);
        x = [0:1/255:1]';
        gamma = [x.^2.2 x.^2.2 x.^2.2];
    end