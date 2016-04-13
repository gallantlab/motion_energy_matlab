function idstr = refmat_store(m,sDir,varNm)
% Usage: idstr = refmat_store(m,sDir,varNm)
% 
% Store matrix temporarily in (sDir) folder of hard drive
% 
% SN 200X?

if ~exist('sDir','var')
    sDir = '/tmp/';
end
if ~exist('varNm','var')
    varNm = 'm';
end
ms = numel(m);
x1 = m(1:ceil(ms/127):end);
tmp.x1 = x1;
tmp.now = now;
idstr = getIDstring(tmp);
mio = matfile(fullfile(sDir,sprintf('refmat_%s.mat', idstr)),'writable',true); 
mio.(varNm) = m;
