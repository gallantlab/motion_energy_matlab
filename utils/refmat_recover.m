function m = refmat_recover(idstr,sDir)
% Usage: m = refmat_recover(idstr,sDir)
% 
% Load a matrix back into working memory after temporary storage on hard
% drive /tmp/ folder
% 
% SN 200X?
if ~exist('sDir','var')
    sDir = '/tmp/';
end
load(fullfile(sDir,sprintf('refmat_%s.mat', idstr)));
