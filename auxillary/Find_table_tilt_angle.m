function [tilt_angle, rotation_angle] = Find_table_tilt_angle(Bdata, voxelSize, limits,offset,showFlag,Nshow)
% Finds and returns the table tilt angle from a set of threshold images.

if nargin < 6
    Nshow = 20;
end

if nargin < 5
    showFlag = true;
end
if isempty(showFlag)
    showFlag = true;
end

if nargin < 4
    offset = 10;
end
if isempty(offset)
    offset = 10;
end

dims = size(Bdata);

if nargin < 3
    limits = [1 dims(3)];
end
if isempty(limits)
    limits = [1 dims(3)];
end

yC = zeros(1,double(dims(3)));
rotation_angles = yC;

for i = 1:dims(3)
    % build L/R masks
    BdataL = 0*Bdata(:,:,i);
    BdataR = BdataL;
    BdataL(:,1:floor(dims(2)/2)-offset) = Bdata(:,1:floor(dims(2)/2)-offset,i);
    BdataR(:,ceil(dims(2)/2)+offset:dims(2)) = Bdata(:,ceil(dims(2)/2)+offset:dims(2),i);
    BdataL = (BdataL>0.5);    
    BdataR = (BdataR>0.5);
    
    % compute coordinates of lowest point
    sumBdataL = sum(BdataL,2);
    yL = dims(1);
    stopFlag = false;
    for j = (dims(1)-1):-1:1
        if sumBdataL(j+1) == 0 && sumBdataL(j) > 0 && ~stopFlag
            yL = j;
            stopFlag = true;
        end
    end
    xLs = double(find(BdataL(yL,:)));
    xL = sum(xLs)/length(xLs);
            
    sumBdataR = sum(BdataR,2);
    yR = dims(1);
    stopFlag = false;
    for j = (dims(1)-1):-1:1
        if sumBdataR(j+1) == 0 && sumBdataR(j) > 0 && ~stopFlag
            yR = j;
            stopFlag = true;
        end
    end
    xRs = double(find(BdataR(yR,:)));
    xR = sum(xRs)/length(xRs);

    m = (yR - yL)/(xR - xL);
    rotation_angles(i) = atan2d((yR - yL),(xR - xL));
    xC = (double(dims(2))+1)/2.0;
    yC(i) = m*(xC-xL) + yL;    
    
    if mod(i,Nshow) == 0 && showFlag
        figure, 
        imshow(Bdata(:,:,i),[]);
        line([xL,xC],[yL,yC(i)],...
               'LineWidth',3,...
               'Color',[.2,.5,1]); 
        line([xC,xR],[yC(i),yR],...
               'LineWidth',3,...
               'Color',[.2,.5,1]);  
        drawnow;
    end        
end
slices = 1:dims(3);


if showFlag
figure;
plot(slices,yC);
xlabel('Scan slice number')
ylabel('y-coord of center of table')
title('Y-center vs slice')

figure;
plot(slices,rotation_angles);
xlabel('Scan slice number')
ylabel('Rotation angles (degrees)')
title('Rotation angle vs slice')
end

z = slices(limits(1):limits(2))';
y = yC(limits(1):limits(2))';

% perform linear regression y = beta(1) + beta(2)*z 
Z = [ones(length(z),1),z];
beta = Z\y; 

rise = voxelSize(1)*beta(2)*(limits(2)+1-limits(1));
run = (limits(2)+1-limits(1))*voxelSize(3);

tilt_angle = atan2d(rise,run);

% extremum correction
ymax = 0;
if tilt_angle > 0
    for k = z(1):z(end)
        if yC(k) < ymax
            y(k-z(1)+1) = ymax;
        else
            ymax = yC(k);
        end
    end
else
    for k = z(end):-1:z(1)
        if yC(k) < ymax
            y(k-z(1)+1) = ymax;
        else
            ymax = yC(k);
        end
    end    
end

% perform linear regression y = beta(1) + beta(2)*z 
beta = Z\y; 
yFit = Z*beta;

rise = voxelSize(1)*beta(2)*(limits(2)+1-limits(1));
run = (limits(2)+1-limits(1))*voxelSize(3);

tilt_angle = atan2d(rise,run);

figure,
scatter(z,y)
hold on
plot(z,yFit)
daspect([voxelSize(1) voxelSize(3) voxelSize(1)])
xlabel('Scan slice number')
ylabel('Table coordinate')
title('Table coordinate vs slice')
grid on

%size(rotation_angles)
rotation_angle = nanmean(rotation_angles,2);

end

