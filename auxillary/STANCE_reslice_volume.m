function [V_new, Y_new] = STANCE_reslice_volume(V_old,scan,sumThreshold,ZcenterFlag,shift2GridFlag,method)
%%% Spontaneous and Task-related Activation of Neurons with Correlations Emulator (STANCE) %%%
%STANCE__reslice_volume generates a resliced volume of V_old
%   using scan parameters into V_new:
%   sumThreshold is an integer that sets the requirement for a slice to be
%   included; if an array it specifies the slice limits
%
%        V = same format as SPM12 header
%        Y = MRI data array
%
%        scan.voxel.size: 3x1 list of voxel physical dimensions [mm]
%                        (default) = [3 3 3]
%        scan.voxel.dimensions: 3x1 list of the number of voxels per spatial direction
%                        (default) = [64 64 NaN]
%                        NOTE: dimensions(3) = NaN used indicate auto
%                              determine setting
%        scan.voxel.spacing: 3x1 list of inter-voxel spacing [mm]
%        scan.tiltAngle: tilt angle in YZ plane away from nose [degrees]
%                        (default) = 15
%
%        sumThreshold: (option 1) numerical threshold required upon a slice 
%                                 for inclusion in the resulting 3D volume.
%                      (option 2) a 2 element array where
%                                 sumThreshold(1) = the lower slice index
%                                 to use after the affine transformation
%                                 sumThreshold(1) = the upper slice index
%                                 to use after the affine transformation   
%
%        ZcenterFlag:   true then use the volume center not the origin
%                       (default) = false
%
%        method: sets the interpolation method for the resampling.
%                       0          Zero-order hold (nearest neighbour).
%            (default)  1          First-order hold (trilinear interpolation).
%                       2->127     Higher order Lagrange (polynomial) interpolation using
%                                  different holds (second-order upwards).
%                      -127 - -1   Different orders of sinc interpolation.
%
% Jason E. Hill
% STANCE_reslice_volume.m      updated     14 MAR 2016

%% Input argument handling and initialization

if nargin < 6
    method = 1;
end
if nargin < 5
    shift2GridFlag = false;
end
if nargin < 4
    ZcenterFlag = false; 
end
if nargin < 3
    sumThreshold   = 100;
    sliceLimitFlag = false; %#ok<NASGU>
end
if length(sumThreshold) == 2
    sliceLimitLower = sumThreshold(1);
    sliceLimitUpper = sumThreshold(2);
    sliceLimitFlag  = true;
else
    sliceLimitFlag  = false;
end    

if nargin < 2
    scan = [];
end

dimMax = 320; % limit any dimension to this for sake of memory consumption/error handling

%% Load scan info 

if nargin < 2 || isempty(scan)
    voxelSize    = [3 3 3];
    new_dims     = [64 64 NaN];  % [64 64 40];
    tiltAngle    = 15; % degrees
    voxelSpacing = [0 0 0.6];
else
    voxelSize    = scan.voxel.size;
    new_dims     = scan.voxel.dimensions;
    voxelSpacing = scan.voxel.spacing;    
    tiltAngle    = scan.tiltAngle;
end
if isnan(new_dims(3))
    sliceNbrFlag = false;
else
    sliceNbr = new_dims(3);
    sliceNbrFlag = true;
end

dX = voxelSize(1) + voxelSpacing(1);
dY = voxelSize(2) + voxelSpacing(2);
dZ = voxelSize(3) + voxelSpacing(3);

%% Load header info

old_dims = [V_old.dim(1), V_old.dim(2), V_old.dim(3)]';
originH = V_old.mat\[0 0 0 1]';
origin = originH(1:3);
dX_old = abs(V_old.mat(1,1));
dY_old = abs(V_old.mat(2,2));
dZ_old = abs(V_old.mat(3,3));
%origin = [origin(1)/dX_old,origin(2)/dY_old,origin(3)/dZ_old]';  
    
% convert tilt angle to "grid" value
Delta_Y = (dZ_old/dY_old)*tand(tiltAngle);
tiltAngle = atan2d(Delta_Y,1);

%old_dims
origin
if sum(old_dims-origin) < double(sum(old_dims))/30.0;
    error('The origin is not near the center of the image!')
end
    
%% determine necessary values for affine transformation

Xc = (old_dims(1) + 1)/2;
Yc = (old_dims(2) + 1)/2;

if mod(old_dims(3),2)
    % odd case
    Zc = (old_dims(3) + 1)/2;
    evenFlag = false;
else
    Zc = (old_dims(3))/2;
    evenFlag = true;
end

Ct = cosd(tiltAngle);
St = sind(tiltAngle);

if shift2GridFlag
  Zout = (origin(3) - Zc)*Ct + (Yc - origin(2))*St;
  if (ceil(Zout)-(Zout))>((Zout) - floor(Zout))
      Zout = floor(Zout);
  else
      Zout = ceil(Zout);
  end
  Zc = origin(3) - (Zout - (Yc - origin(2))*St)/Ct;
end

nZ = ceil((Zc*Ct + Yc*abs(St))*(dZ_old/dZ));

if ~ZcenterFlag
    Zc = origin(3);
end

%% reslice slice by slice according to affine transformation

startFlag  = true;
stopFlag   = false;

if evenFlag
    nMax = 2*nZ;
else
    nMax = (2*nZ+1);
end
nMax = min(nMax,dimMax);

for n = 1:nMax
    if ~stopFlag
        % set rotation about x-axis thru point (x,Yc,Zc)
        % keeping in plane dimensions the same
        %  -> new re-sliced centers along the line: Zk = Zc + T(Yk - Yc),
        %  where the re-sliced center points are at
        %    Ykc = Yc*(1-Ct) + k*St*dZ/dZ_old
        %    Zkc = Zc + k*Ct*dZ/dZ_old
        %     k = -nZ ... 0 ... nZ  for odd # of original slices
        %     k = -nZ+1 ... 0 ... nZ  for even # of original slices        
        %     n = k + (nZ + 1)
        if evenFlag
            k = n - nZ;
        else
            k = n - (nZ+1);
        end
        % re-sliced origin:
        Yk = Yc*(1 - Ct) - k*St*dZ/dZ_old;
        Zk = Zc - Yc*St + k*Ct*dZ/dZ_old; 
        Affine(:,:,n) = [[1 0   0   0]; ...
                         [0 Ct -St  Yk]; ...
                         [0 St  Ct  Zk]; ...                
                         [0 0   0   1]];           %#ok<*AGROW>
        affine_dims = [old_dims(1), old_dims(2)];
        AffineSlice = spm_slice_vol(V_old,Affine(:,:,n),affine_dims,method);    
        
        slicesSums(n) = squeeze(sum(sum(AffineSlice)));
%        figure, imshow(AffineSlice);
        if  ~sliceLimitFlag
            if sliceNbrFlag
                AffineSlices(:,:,n) = AffineSlice;
            elseif (slicesSums(n) > sumThreshold) && startFlag && ~stopFlag
                startFlag = false;
                AffineSlices(:,:,1) = AffineSlice;
                slice_n = 1;
                sliceLimitLower = n;
            elseif (slicesSums(n) > sumThreshold) && ~startFlag && ~stopFlag
                slice_n = slice_n + 1;
                AffineSlices(:,:,slice_n) = AffineSlice;
%                figure, imshow(sA);
            elseif (slicesSums(n) < sumThreshold) && ~startFlag && ~stopFlag
                stopFlag = true;
                sliceLimitUpper = n-1;
            end
        else
            if n == sliceLimitLower
                slice_n = 1;
                AffineSlices(:,:,1) = AffineSlice;
                startFlag = false;
            elseif (n > sliceLimitLower)&&(n < sliceLimitUpper)
                slice_n = slice_n + 1;
                AffineSlices(:,:,slice_n) = AffineSlice;
            end
            if n == sliceLimitUpper
                slice_n = slice_n + 1;
                AffineSlices(:,:,slice_n) = AffineSlice;
                stopFlag = true;
            end 
        end
    end
end
if sliceNbrFlag    
%    sliceNbr
    [~, idxSort]  = sort(slicesSums,'descend');
%    slicesSums
%    idxSort
    sliceLimit = idxSort(sliceNbr);
    if sliceLimit < (nZ + 1)
        sliceLimitLower = sliceLimit;
        if sliceLimit + sliceNbr - 1 > nMax
            sliceLimitUpper = nMax;
            sliceLimitLower = sliceLimit - (nMax - (sliceLimit + sliceNbr - 1));
            if sliceLimitLower < 1
                sliceLimitLower = 1;
                sliceNbrWarning = sprintf('Cannot accomodate the requested number of slices, using %d slices.',nZ);
                warning(sliceNbrWarning); %#ok<*SPWRN>
            end
        else
            sliceLimitUpper = sliceLimit + sliceNbr - 1;
        end
    else
        sliceLimitUpper = sliceLimit;
        if sliceLimit - sliceNbr + 1 < 1
            sliceLimitLower = 1;
            if sliceNbr > nMax     
                slice = nMax; %#ok<NASGU>
                sliceNbrWarning = sprintf('Cannot accomodate the requested number of slices, using %d slices.',nZ);
                warning(sliceNbrWarning);                
            else
                sliceLimitUpper = sliceNbr;
            end
        else
            sliceLimitLower = sliceLimit - sliceNbr + 1;
        end
    end
    AffineSlices = AffineSlices(:,:,sliceLimitLower:sliceLimitUpper);
    slice_n = size(AffineSlices,3);
end
affine_dims = size(AffineSlices);
affine_dims(1) = min(affine_dims(1),dimMax);
affine_dims(2) = min(affine_dims(2),dimMax);

AffineOrigin(1) = origin(1);
AffineOrigin(2) = Yc + (origin(2)-Yc)*Ct + (origin(3)-Zc)*St;
if evenFlag
    AffineOrigin(3) = (nZ - sliceLimitLower + 1) + (origin(3) - Zc)*Ct*(dZ_old/dZ) + (Yc - origin(2))*St*(dZ_old/dZ);
else
    AffineOrigin(3) = (nZ - sliceLimitLower + 2) + (origin(3) - Zc)*Ct*(dZ_old/dZ) + (Yc - origin(2))*St*(dZ_old/dZ);    
end
%AffineOrigin

%% reslice to target scan dimensions

newImageCenter(1) = 0.5*(affine_dims(1)+1);
newImageCenter(2) = 0.5*(affine_dims(2)+1);

% if ~ZcenterFlag
newOrigin(1) = newImageCenter(1)*dX + (AffineOrigin(1) - Xc)*dX_old;
newOrigin(2) = newImageCenter(2)*dY + (AffineOrigin(2) - Yc)*dY_old;
newOrigin(3) = AffineOrigin(3)*dZ;
% else
% newOrigin(1) = newImageCenter(1)*dX + (AffineOrigin(1) - Xc)*dX_old;
% newOrigin(2) = newImageCenter(2)*dY + (AffineOrigin(2) - Yc)*dY_old;
% newOrigin(3) = AffineOrigin(3)*dZ;    
% end

newOrigin

[YA,XA] = meshgrid(1:affine_dims(2),1:affine_dims(1));
[YN,XN] = meshgrid(newImageCenter(2)-0.5*(dY/dY_old)*(new_dims(2)-1):(dY/dY_old):newImageCenter(2)+0.5*(dY/dY_old)*(new_dims(2)-1),...
                   newImageCenter(1)-0.5*(dX/dX_old)*(new_dims(1)-1):(dX/dX_old):newImageCenter(1)+0.5*(dX/dX_old)*(new_dims(1)-1));
               
newData = zeros(new_dims(1),new_dims(2),slice_n);
%size(newData)
for n = 1:slice_n  
    newData(:,:,n) = interp2(YA,XA,AffineSlices(:,:,n),YN,XN);
%    figure, imshow(newData(:,:,n));
end

%% save to new volume structure
fileName = V_old.fname;
fileExtension = fileName(end-2:end);
if strcmpi(fileExtension,'.gz')
    fileExtension = fileName(end-6:end-3);
    fileBase = fileName(1:end-7);
else
    fileBase = fileName(1:end-4);
    fileExtension = fileName(end-3:end);
end
dataFileNameNIIout = [fileBase,'_resliced',fileExtension];

V_new.fname = dataFileNameNIIout;
V_new.dim = size(newData);
V_new.dt = V_old.dt;
V_new.pinfo = V_old.pinfo;
V_new.n = [1,1];
V_new.descrip = V_old.descrip;
V_new.private = V_old.private;
Y_new = newData;
Y_new(isnan(Y_new)) = 0;

if length(sumThreshold) == 1
    V_new.sliceLimitLower = sliceLimitLower;
    V_new.sliceLimitUpper = sliceLimitUpper;
end

% generate new world coordinate matrix
%V_old.mat

sign_mat = sign(V_old.mat);
new_mat = zeros(4);
new_mat(1,1) = sign_mat(1,1)*dX;
new_mat(2,2) = sign_mat(2,2)*dY;
new_mat(3,3) = sign_mat(3,3)*dZ;
new_mat(4,4) = 1; 

new_mat(1,4) = sign_mat(1,4)*abs(newOrigin(1));
new_mat(2,4) = sign_mat(2,4)*abs(newOrigin(2));
new_mat(3,4) = sign_mat(3,4)*abs(newOrigin(3));

%new_mat

V_new.mat = new_mat;

V_new = spm_write_vol(V_new,Y_new);

end

