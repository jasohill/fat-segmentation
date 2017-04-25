function [x_cm, y_cm, angle_xz, angle_yz] = Find_cm_angle(Bdata, voxelSize, slice_thickness,limits,offset,showFlag)
% Finds and returns the angle of the plane relative to the center line of the table from a set of threshold images.

if nargin < 6
    showFlag = true;
end
if isempty(showFlag)
    showFlag = true;
end

if nargin < 5
    offset = 10;
end
if isempty(offset)
    offset = 10;
end

dims = size(Bdata);

if nargin < 4
    limits = [1 dims(3)];
end
if isempty(limits)
    limits = [1 dims(3)];
end

x_cm = zeros(1,double(dims(3)));
y_cm = zeros(1,double(dims(3)));


for k = 1:dims(3)
    % find center of mass
    CoM = center_of_mass(Bdata(:,:,k)>0.5);
    y_cm(k) = CoM(1);
    x_cm(k) = CoM(2);
end
slices = 1:dims(3);

if showFlag
figure;
plot(slices,x_cm);
title('Center of mass (x-coord) vs. slice number')

% figure;
% plot(slices,y_cm);
% title('Center of mass (y-coord) vs. slice number')
end

% find the angle that the body makes in the XZ plane
z = slices(limits(1):limits(2))';
x = x_cm(limits(1):limits(2))';

% perform linear regression x = beta(1) + beta(2)*z 
Z = [ones(length(z),1),z];
beta = Z\x; 
xFit = Z*beta;

rise = voxelSize(1)*beta(2)*(limits(2)+1-limits(1));
run = (limits(2)+1-limits(1))*slice_thickness;

angle_xz = atan2d(rise,run);

if showFlag
figure,
scatter(z,x)
hold on
plot(z,xFit)
daspect([voxelSize(2) slice_thickness voxelSize(1)])
xlabel('Scan slice number')
ylabel('x coordinate of CM')
title('x coordinate of CM vs slice')
grid on

% find the angle that the body makes in the YZ plane
z = slices(limits(1):limits(2))';
y = y_cm(limits(1):limits(2))';

% perform linear regression x = beta(1) + beta(2)*z 
Z = [ones(length(z),1),z];
beta = Z\y; 
yFit = Z*beta;

rise = voxelSize(1)*beta(2)*(limits(2)+1-limits(1));
run = (limits(2)+1-limits(1))*slice_thickness;

angle_yz = atan2d(rise,run);

% if showFlag
%  figure,
%  scatter(z,y)
%  hold on
%  plot(z,yFit)
%  daspect([voxelSize(2) slice_thickness voxelSize(1)])
%  xlabel('Scan slice number')
%  ylabel('y coordinate of CM')
%  title('y coordinate of CM vs slice')
%  grid on
%  end

end