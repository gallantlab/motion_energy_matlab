% Demo for colorspace.m - the CIE xyY "tongue"

% Pascal Getreuer 2006

N = 160;
Nx = round(N*0.8);
Ny = round(N*0.9);
e = 0.01;
% Generate colors in the xyY color space
x = linspace(e,0.8-e,Nx);
y = linspace(e,0.9-e,Ny);
[xx,yy] = meshgrid(x,y);
iyy = 1./(yy + 1e-5*(yy == 0));

% Convert from xyY to XYZ
Y = ones(Ny,Nx);
X = Y.*xx.*iyy;
Z = Y.*(1-xx-yy).*iyy;
% Convert from XYZ to R'G'B'
CIEColor = colorspace('rgb<-xyz',cat(3,X,Y,Z));

figure;
% Render the colors on the tonue
image(x,y,CIEColor)
colormap(gray);
set(gca,'YDir','normal');
axis image

% The boundary of the tongue
xy = [[0.1740, 0.0050];[0.1736, 0.0049];[0.1703, 0.0058];
   [0.1566, 0.0177];[0.1440, 0.0297];[0.1241, 0.0578];
   [0.1096, 0.0868];[0.0913, 0.1327];[0.0687, 0.2007];
   [0.0454, 0.2950];[0.0235, 0.4127];[0.0082, 0.5384];
   [0.0039, 0.6548];[0.0139, 0.7502];[0.0389, 0.8120];
   [0.0743, 0.8338];[0.1142, 0.8262];[0.1547, 0.8059];
   [0.1929, 0.7816];[0.2296, 0.7543];[0.2658, 0.7243];
   [0.3016, 0.6923];[0.3373, 0.6589];[0.3731, 0.6245];
   [0.4087, 0.5896];[0.4441, 0.5547];[0.4788, 0.5202];
   [0.5125, 0.4866];[0.5448, 0.4544];[0.5752, 0.4242];
   [0.6029, 0.3965];[0.6270, 0.3725];[0.6482, 0.3514];
   [0.6658, 0.3340];[0.6801, 0.3197];[0.7006, 0.2993];
   [0.7140, 0.2859];[0.7260, 0.2740];[0.7340, 0.2660]];
% Make a smooth boundary with spline interpolation
xi = [interp1(xy(:,1),1:0.25:size(xy,1),'spline'),xy(1,1)];
yi = [interp1(xy(:,2),1:0.25:size(xy,1),'spline'),xy(1,2)];

% Draw the boundary of the tongue
hold on;
set(patch([0.8,0,0,0.8,xi(:).',0.8],[0.9,0.9,0,0,yi(:).',0],...
   [1,1,1]),'EdgeColor','none');
plot(xi,yi,'k-');

axis([0,0.8,0,0.9]);
xlabel('x');
ylabel('y');
title('The CIE "tongue": the region of all colors over x and y');
shg;
