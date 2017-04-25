function [x_cm, angleL, angleR] = Find_plane_angle(Bdata, voxelSize, slice_thickness,limits,offset,showFlag)
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

xL = zeros(1,double(dims(3)));
xR = zeros(1,double(dims(3)));

%yC = zeros(1,double(dims(3)));
%rotation_angles = yC;

for k = 1:dims(3)
    % build L/R masks
    BdataL = 0*Bdata(:,:,k);
    BdataR = BdataL;
    BdataL(:,1:floor(dims(2)/2)-offset) = Bdata(:,1:floor(dims(2)/2)-offset,k);
    BdataR(:,ceil(dims(2)/2)+offset:dims(2)) = Bdata(:,ceil(dims(2)/2)+offset:dims(2),k);
    BdataL = (BdataL>0.5);    
    BdataR = (BdataR>0.5);

    % compute coordinates of leftmost point
    sumBdataL = sum(BdataL,1);
    xL(k) = 1;
    stopFlag = false;
    for i = 2:dims(2)
        if sumBdataL(i-1) == 0 && sumBdataL(i) > 0 && ~stopFlag
            xL(k) = i;
            stopFlag = true;
        elseif sumBdataL(i-1) > 0
            stopFlag = true;
        end
    end
    yLs = double(find(BdataL(:,xL(k))));
    yL(k) = sum(yLs)/length(yLs);    

    % compute coordinates of rightmost point
    sumBdataR = sum(BdataR,1);
    xR(k) = dims(2);
    stopFlag = false;
    for i = (dims(2)-1):-1:1
        if sumBdataL(i+1) == 0 && sumBdataR(i) > 0 && ~stopFlag
            xR(k) = i;
            stopFlag = true;
        elseif sumBdataL(i+1) > 0
            stopFlag = true;
        end
    end
    yRs = double(find(BdataL(:,xR(k))));
    yR(k) = sum(yRs)/length(yRs);                 
end
slices = 1:dims(3);

if showFlag
figure;
plot(slices,xL);
title('Left side of body vs. slice number')

figure;
plot(slices,xR);
title('Right side of body vs. slice number')
end

% find the x center of mass
x_cm = 0.5*(mean(xL) + mean(xR));

% find the angle that the body makes on the left
z = slices(limits(1):limits(2))';
x = xL(limits(1):limits(2))';

% perform linear regression x = beta(1) + beta(2)*z 
Z = [ones(length(z),1),z];
beta = Z\x; 
xFit = Z*beta;

rise = voxelSize(1)*beta(2)*(limits(2)+1-limits(1));
run = (limits(2)+1-limits(1))*slice_thickness;

angleL = atan2d(rise,run);

if showFlag
figure,
scatter(z,x)
hold on
plot(z,xFit)
daspect([voxelSize(2) slice_thickness voxelSize(1)])
xlabel('Scan slice number')
ylabel('x coordinate of left side')
title('x coordinate left vs slice')
grid on
end

% find the angle that the body makes on the right
x = xR(limits(1):limits(2))';

% perform linear regression x = beta(1) + beta(2)*z 
Z = [ones(length(z),1),z];
beta = Z\x; 
xFit = Z*beta;

rise = voxelSize(1)*beta(2)*(limits(2)+1-limits(1));
run = (limits(2)+1-limits(1))*slice_thickness;

angleR = atan2d(rise,run);

if showFlag
figure,
scatter(z,x)
hold on
plot(z,xFit)
daspect([voxelSize(2) slice_thickness voxelSize(1)])
xlabel('Scan slice number')
ylabel('x coordinate of right side')
title('x coordinate right vs slice')
grid on
end

end

