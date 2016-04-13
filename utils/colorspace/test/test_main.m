% Test script for colorspace.m

% Pascal Getreuer 2006

%%% Conversion Accuracy Test %%%

fprintf(['\n Conversion Accuracy Test\n\n',...
      'To verify the invertibility of the color transfomations, this test\n',...
      'converts R''G''B'' data to a space, inverts, and compares with the\n',...
      'original data.\n']);
N = 1e4;            % Number of points to test
A = rand(N,1,3);    % Generate points in R'G'B' colorspace

Space = {'YPbPr','YCbCr','YDbDr','JPEG-YCbCr','YIQ','YUV',...
      'HSV','HSL','XYZ','Lab','Luv','Lch'};
fprintf('\n Conversion        Max Error\n\n');

for k = 1:length(Space)
   B = colorspace([Space{k},'<-RGB'],A);  % Convert to Space{k}
   R = colorspace(['RGB<-',Space{k}],B);  % Convert back to R'G'B'
   MaxError = max(abs(A(:) - R(:)));
   fprintf(' RGB<->%-10s  %.2e\n',Space{k},MaxError);
end

fprintf('\n\n');