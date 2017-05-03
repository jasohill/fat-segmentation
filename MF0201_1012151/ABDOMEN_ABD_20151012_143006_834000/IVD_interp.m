% Authors:     Ryan Hayden & Jason E. Hill, Ph.D.,
% Institution: Texas Tech University
%              Dept. of Computer and Electrical Engineerin
% Updated:     19 APR 2017

clearvars, close all;
currDir = pwd;
cd ..
cd ..
addpath(genpath(pwd))
cd(currDir)

%% load volumes
load subject01_i_registered.mat
resultsFile = 'subject01_i_IVDinterp.mat';
writeResults = true;
showProgress = false;
showSlice = 53;

niiWup_reslice = load_nii(NIFTI_file_name_W_upper_reslice);
niiWlo_reslice = load_nii(NIFTI_file_name_W_lower_reslice);

% slice loop prep
indeces = 1:Nslices;
%first_i  = indeces(1);
DW = 320; % Slice width
DH = 260; % Slice height

%% Extract water image data; inhomogeneity correct; auto-threshold 
h = waitbar(0,'Please wait, enhancing water image data...');
WaterImages = uint8(zeros(DH,DW,Nslices));
for i = indeces
   slice = slices(i);   
   level = levels(i); 

   % make 8-bit slice image of water selected data
   % orient image slice (posterior down)
   if level == 1
       IWraw = fliplr(rot90(niiWup_reslice.img(:,:,slice)));    
   elseif level == 2
       IWraw = fliplr(rot90(niiWlo_reslice.img(:,:,slice)));       
   end

   IWrawSample = 0*IWraw;
   IWrawSample(50:200,90:230) = IWraw(50:200,90:230);
   IWrawSample(110:135,135:175) = mean(mean(IWraw));
   IWrawSampleMax = max(max(IWrawSample));
   IWraw(IWraw > IWrawSampleMax) = IWrawSampleMax;
   % convert data to 8-bit gray-scale image
   IWshow = uint8(double(IWraw)/double(IWrawSampleMax)*255);

   IWrawVeto = (IWshow < vetoFactorW*backgroundThreshold);
   IWrawVeto = logical((1 - ~(IWshow < 0.75*vetoFactorW*backgroundThreshold).*imdilate(~IWrawVeto,ones(5))));
   IWrawVeto = logical((1 - ~(IWshow < backgroundThreshold).*imdilate(~IWrawVeto,ones(3))));
        
   %% correct water signal inhomogeneity

   % construct rough estimate of the signal inhomogeneity via morphological closing
   [IWshow2,~] = CorrectInhomogeneity(IWraw,50,[],1.0); %50 %100);

if showProgress && i == showSlice    
figure,
imshow(IWshow2)   
end
   
   %% segment the image foreground from background within body cavity 

   % rescale image to the maximum intensity
   IWmax = max(max(IWshow2));
   IWshow3 = uint8(double(IWshow2)/double(IWmax)*255);
   IWshow3 = IWshow3.*uint8(1-IWrawVeto);


if showProgress && i == showSlice   
figure,
imshow(IWshow3)   
end
   
   % Find the peak and valley pixel bins of histogram
if showProgress && i == showSlice  
   [peakPixelBinW,valleyPixelBinW] = FindHistogramExtrema(IWshow3,backgroundThreshold,foregroundThresholdW+43,true);
   peakPixelBinW = peakPixelBinW + 2*backgroundThreshold %#ok<NOPTS>
   peakIndicesW = find(IWshow3 == peakPixelBinW,1) %#ok<NOPTS>
else
   [peakPixelBinW,valleyPixelBinW] = FindHistogramExtrema(IWshow3,backgroundThreshold,foregroundThresholdW+43,false);
end
   peakPixelBinW = peakPixelBinW + 2*backgroundThreshold;
   if peakPixelBinW > 255
       peakPixelBinW = 255;
   end
   %valleyPixelBinW = valleyPixelBinW %- backgroundThreshold

   % rescale water image using the histogram peak as the maximum  
   peakIndicesW = find(IWshow3 == peakPixelBinW);
   if isempty(peakIndicesW)
       peakIndicesW = find(IWshow3 == peakPixelBinW + 1);
   end

   IWpeak = IWshow2(peakIndicesW(1));
   IWshow2(IWshow2>IWpeak) = IWpeak;
   IWshow2 = IWshow2.*uint8(1-IWrawVeto);

if showProgress && i == showSlice   
figure,
imshow(IWshow2)   
end  
   
   WaterImage = uint8(double(IWshow2)/double(IWpeak)*255);

   % limit the valley pixel bin by Otsu's multithreshold method
   valleyPixelBinW = AutorestrictThreshold(WaterImage,valleyPixelBinW,nbrThresholdsW);

   % define water foreground as being above histogram valley
   BinaryWaterImage = logical((WaterImage > valleyPixelBinW).*double(1-IWrawVeto));

   WaterImages(:,:,i) = rot90(WaterImage,2).*uint8(rot90(BinaryWaterImage,2));
   
   waitbar(i/Nslices)
   
end
close(h)

% RESULT: a WaterImage volume & a BinaryWaterImage

%% Interpolation Grid
IVD.zs     = [75, 90, 103, 115, 10, 23, 34, 46, 58, 67, 76, 84, 95, 102, 110];
IVD.levels = [2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
IVD.labels = [-4, -3, -2, -1, -1, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3];
IVD.xs     = [164, 164, 164, 164, 161, 163, 163, 164, 164, 164, 164, 167, 167,  167,  167];  
IVD.ys     = [137, 139, 133, 128, 129, 124, 119, 112, 109, 106, 106, 105, 104,  105,  106];      

IVD.N = length(IVD.zs);

IVD.zs_interp = zeros(1,IVD.N);
for k = 1:IVD.N
    
z = IVD.zs(k);

% apply to slice
slice = get_slice_index(slices,levels,z,IVD.levels(k)); %#ok<NOPTS>

WI = WaterImages(:,:,slice);

% initial IVD mask seed
W = 11; % Width of intervertebral disc mask
H = 7;  % Height of intervertebral disc mask 

% Middle of intevertebral disc coordinates
x = IVD.xs(k); 
y = IVD.ys(k); 

stop_flag = false;
stop_threshold = 50;
max_iteration  = 25;
iter = 1;
while ~stop_flag && (iter < max_iteration)
    [X, Y] = meshgrid(1:DW,1:DH); % Image size 
    EMask = ((X - x)/W).^2 + ((Y - y)/H).^2 <= 1;
    BMask = (EMask.*logical(WI>0));
    mask_sum = sum(EMask(:)) - sum(BMask(:));
    if mask_sum < stop_threshold
        % expand mask
        W = W+1;
        H = H+1;
        % update the center of mass
        IVDmask = WI.*uint8(EMask);
        CoM = center_of_mass(IVDmask);
        X = CoM(1);
        Y = CoM(2);
    else
        stop_flag = true;
    end
    iter = iter+1;
end

% clean up mask with morphological operations
se3 = strel('diamond',1);
BMask = imopen(BMask,se3);
% se11 = strel('disk',5);
BMask = imopen(imfill(BMask),se3);
BMask = imclose(BMask,se3);

%% Apply opened ellipse mask to slice to creat the IVD mask
figure,
imshow(WI)

IVDmask = WI.*uint8(BMask);
figure,
imshow(IVDmask)

% Determine the intensities of the neighboring slices
stop_flag = false;
while ~stop_flag
   
   if slice > 2 && levels(slice-1) == levels(slice)
       WI_minus_1 = WaterImages(:,:,slice-1);
       IVDmask_minus_1 = WI_minus_1.*uint8(EMask);
   else
       IVDmask_minus_1 = IVDmask;
   end
   if slice < (Nslices-1) && levels(slice+1) == levels(slice)
       WI_plus_1 = WaterImages(:,:,slice+1);
       IVDmask_plus_1 = WI_plus_1.*uint8(EMask);       
   else
       IVDmask_plus_1 = IVDmask;
   end
   if sum(IVDmask_minus_1(:)) > sum(IVDmask(:))
       slice = slice - 1;
       z = z + 1;
       IVDmask = IVDmask_minus_1;
   elseif sum(IVDmask_plus_1(:)) > sum(IVDmask(:))
       slice = slice + 1;
       z = z - 1;
       IVDmask = IVDmask_plus_1;       
   else
       stop_flag = true;
   end
end


deltaZ = [z-1 z z+1];
I = [sum(IVDmask_plus_1(:)) sum(IVDmask(:)) sum(IVDmask_minus_1(:))];
p = polyfit(deltaZ, I, 2);
IVD_z_interp = -p(2)/(2*p(1));
disp(['The interpolated z slice position of the IVD is ',num2str(IVD_z_interp),'.'])

IVD.zs_interp(k) = IVD_z_interp; 

IVD.mask(:,:,k) = IVDmask;

if isempty(get_slice_index(slices,levels,round(IVD_z_interp),IVD.levels(k)))
   IVD.slice(k) = get_slice_index(slices,levels,z,IVD.levels(k));
else
   IVD.slice(k) = get_slice_index(slices,levels,round(IVD_z_interp),IVD.levels(k));
end

end

IVD.zs_interp

% correct scan offset
IVD.scanOffset = IVD.zs_interp(4)-IVD.zs_interp(5);   %Z_L1_L2 - Z_L1_L2u (98 discrete)

IVD.first_lower_l = 121 - lowerTop + upperBottom + IVD.scanOffset;

%NOTE: add this to a IVD.zs array for all IVDs and 
% also construct a IVD.labels with this one labeled 12 for T12.
% The next one down will be labelled -1 for L1 and so forth.

if writeResults       
    save(resultsFile,'IVD','WaterImages');
end

%% create spinal cord mask
SpineVol = WaterImages;
SpineVol(WaterImages == 0) = 92;
k = 1;
for s = length(indeces):-1:1
   if s < IVD.slice(k)
       if k < length(IVD.slice)
          k = k + 1;
       end
   end
   if k == 1
       SpineVol(:,:,s) =  (SpineVol(:,:,s)).*uint8((IVD.mask(:,:,1)>0));
   else
       SpineVol(:,:,s) =  (SpineVol(:,:,s)).*uint8((IVD.mask(:,:,k-1)>0)|(IVD.mask(:,:,k)>0));
   end
       if ~logical(mod(s,10)) && k > 1
          figure, imshow(squeeze(SpineVol(:,:,s)),[]);
       end   
end
    
%% visualize the spinal cord
scanOffset = round(IVD.scanOffset);

first_lower_i = length(upperSlices) + lowerTop - upperBottom - scanOffset;

show_i = [1:length(upperSlices) first_lower_i:Nslices];
overlap_lower_i = (length(upperSlices)+1):(first_lower_i-1);
overlap_upper_i = (length(upperSlices)-length(overlap_lower_i)+1):length(upperSlices);
no_overlap_i = [1:(length(upperSlices)-length(overlap_lower_i)) first_lower_i:Nslices];

Nshow = length(show_i);
Noverlap = length(overlap_lower_i);

SV = SpineVol;
SV(length(upperSlices)) = SV(length(upperSlices)+1);
for i = 1:size(SV,1)
    for j = 1:size(SV,2)
        x = 1:size(SV,3);
        v = double(squeeze(SV(i,j,:)));
        xq = x + IVD.scanOffset - scanOffset; 
        vq = interp1(x,v,xq);
        SV(i,j,length(upperSlices)+1:end) = vq(length(upperSlices)+1:end);
    end
end
for j = 1:Noverlap
    average_weight = j/(Noverlap+1);
    SV(:,:,overlap_upper_i(j))   = average_weight*SV(:,:,overlap_lower_i(j)) + (1-average_weight)*SV(:,:,overlap_upper_i(j));
end

D = flip(SV(:,:,show_i),3);
Ds = smooth3(D);

figure,
hc = vol3d('cdata',D,'texture','3D');
view(3);  
axis([75 200 50 175 (length(indeces)-IVD.slice(1)-20) inf]);   daspect([1 1 1.40625/3.0])
colormap(bone(256));
alphamap('rampup');