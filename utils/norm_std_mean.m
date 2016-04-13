function [s, s_stds, s_means] = norm_std_mean(s, s_stds, s_means)
% Usage: [s, s_stds, s_means] = norm_std_mean(s, s_stds, s_means)
% 
% Compute std and mean of a sample, z-score using them, return z-scored
% data, stds, and means.
% 
% Essentially the same as matlab's "zscore.m", but with optimization for
% very large matrices.
% 
% Commented by ML 2013.03.07

divnum = 20;
snum = size(s,2);
chunksize = ceil(snum/divnum);

if ~exist('s_means', 'var')
	%% calculate mean
	s_means = mean(s, 1);
end

if ~exist('s_stds', 'var')
	%% calculate std
	if prod(size(s)) < 100000 | size(s,2) < 20
		s_stds = std(s, 0, 1);
	else % to avoid out of memory, calculate stds by chunks
		st = 1;
		for ii=1:divnum
			ed = min([snum st+chunksize-1]);
			s_stds(st:ed) = std(s(:,st:ed), 0, 1);
			st = st + chunksize;
		end
	end
end


%% offset for mean
if prod(size(s)) < 100000 | size(s,2) < 20
    % meanmat = repmat(s_means, [size(s,1) 1]);
    % s = s - meanmat;
    s = bsxfun(@minus, s, s_means);
else
	st = 1;
	for ii=1:divnum
		ed = min([snum st+chunksize-1]);

        % meanmat = repmat(s_means(st:ed), [size(s,1) 1]);
        % s(:,st:ed) = s(:,st:ed) - meanmat;
        s(:,st:ed) =  bsxfun(@minus, s(:,st:ed), s_means(st:ed));

		st = st + chunksize;
	end
end


s_stds(find(s_stds==0)) = 1; % To avoid zero-division

%% normalize stds
if prod(size(s)) < 100000 | size(s,2) < 20
    % stdmat = repmat(s_stds, [size(s,1) 1]);
    % s = s./stdmat;
    s = bsxfun(@rdivide, s, s_stds);
else
	st = 1;
	for ii=1:divnum
		ed = min([snum st+chunksize-1]);
        % stdmat = repmat(s_stds(st:ed), [size(s,1) 1]);
        % s(:,st:ed) = s(:,st:ed)./stdmat;
        s(:,st:ed) = bsxfun(@rdivide, s(:,st:ed), s_stds(st:ed));
		st = st + chunksize;
	end
end
