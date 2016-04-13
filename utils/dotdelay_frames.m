function [chouts choutc] = dotdelay_frames(gs, scw, in)
%function out = dotdelay(kernel, in)
%
% calculate linear responses of a system kernel to an input
%
% INPUT:
% [kernel] = kernel N x D matrix where N is number of channels and
%            D is number of delay lines.
%     [in] = input N x S matrix where S is the number of samples
%
% OUTPUT:
%    [out] = vector of kernel responses to input
%

ktsize = size(scw,2);
ktsize2c = ceil(ktsize/2);
ktsize2f = floor(ktsize/2);
itsize = size(in, 2);

% gouts= (gs'*in)';
% goutc= (gc'*in)';
% 
% souts =  gouts*scw(2,:) + goutc*scw(1,:);
% soutc = -gouts*scw(1,:) + goutc*scw(2,:);

gout = (gs*in)';
outs =  gout(:,1)*scw(2,:) + gout(:,2)*scw(1,:);
outc = -gout(:,1)*scw(1,:) + gout(:,2)*scw(2,:);

% the first half
for ii=1:ktsize2c
    z = zeros(ktsize2c-ii,1);
	outs(:,ii) = [outs(ktsize2c-ii+1:end,ii); z];
	outc(:,ii) = [outc(ktsize2c-ii+1:end,ii); z];
end

% the second half
for ii=1:ktsize2f
	ti = ii+ktsize2c;
    z = zeros(ii,1);
	outs(:,ti) = [z; outs(1:end-ii,ti)];
	outc(:,ti) = [z; outc(1:end-ii,ti)];
end
chouts = sum(outs,2);
choutc = sum(outc,2);
