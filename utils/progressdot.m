function progressdot(ii,s1,s2,t)
% Usage: progressdot(ii,s1,s2,t)
% 
% Keeps tabs on progress through a loop with a display of dots at the
% command line.
% 
% Inputs: 
%   ii - index variable to track 
%   s1 - dot every (s1) (prints dot each (s1) ticks through loop)
%   s2 - new line every (s2) (prints time and new line each (s2) ticks through loop)
%    t - end (max of loop)
% 
% Created by SN, commented by ML 2011.08

persistent t0

if ii==1 || ~exist('t0','var') || isempty(t0)
    t0=clock;
end

if ~exist('t', 'var')
    t=NaN;
end

if mod(ii,s1)==0
    fprintf('.');
end
if mod(ii,s2)==0 || ii==t
    fprintf('%d/%d done. (%s) \n',ii,t, disptime(etime(clock, t0)));
end


function str = disptime(sec)

if sec<60
    str = sprintf('%.1f sec', sec);
elseif sec<3600
    str = sprintf('%.1f min', sec/60);
else
    str = sprintf('%.1f hour', sec/3600);
end

