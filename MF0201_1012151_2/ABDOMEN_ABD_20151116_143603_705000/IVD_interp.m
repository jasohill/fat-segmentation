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
load subject01_f_registered.mat
resultsFile = 'subject01_f_IVDinterp.mat';
writeResults = false;
showSlice = 1;

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
BinaryWaterImages = logical(WaterImages);
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

if i == showSlice   
figure,
imshow(IWshow2)   
end
   
   %% segment the image foreground from background within body cavity 

   % rescale image to the maximum intensity
   IWmax = max(max(IWshow2));
   IWshow3 = uint8(double(IWshow2)/double(IWmax)*255);
   IWshow3 = IWshow3.*uint8(1-IWrawVeto);


if i == showSlice   
figure,
imshow(IWshow3)   
end
   
   % Find the peak and valley pixel bins of histogram
if i == showSlice  
   [peakPixelBinW,valleyPixelBinW] = FindHistogramExtrema(IWshow3,backgroundThreshold,foregroundThresholdW+42,true);
   peakPixelBinW = peakPixelBinW + 2*backgroundThreshold %#ok<NOPTS>
   peakIndicesW = find(IWshow3 == peakPixelBinW,1) %#ok<NOPTS>
else
   [peakPixelBinW,valleyPixelBinW] = FindHistogramExtrema(IWshow3,backgroundThreshold,foregroundThresholdW+45,false);
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

if i == showSlice   
figure,
imshow(IWshow2)   
end  
   
   WaterImage = uint8(double(IWshow2)/double(IWpeak)*255);

   % limit the valley pixel bin by Otsu's multithreshold method
   valleyPixelBinW = AutorestrictThreshold(WaterImage,valleyPixelBinW,nbrThresholdsW);

   % define water foreground as being above histogram valley
   BinaryWaterImage = logical((WaterImage > valleyPixelBinW).*double(1-IWrawVeto));

   WaterImages(:,:,i)       = rot90(WaterImage,2);
   BinaryWaterImages(:,:,i) = rot90(BinaryWaterImage,2);
   %WaterImages(:,:,i) = WaterImage;
   %BinaryWaterImages(:,:,i) = BinaryWaterImage;
   waitbar(i/Nslices)
   
end
close(h)

% RESULT: a WaterImage volume & a BinaryWaterImage

%% Interpolation Grid

IVD.zs     = [69, 83, 98, 110, 11, 23, 35, 47, 58, 68, 77, 86, 95, 103, 111];
IVD.levels = [2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
IVD.labels = [-4, -3, -2, -1, -1, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3];
IVD.xs     = [165, 165, 164, 164, 164, 162, 163, 163, 164, 164, 164, 165, 165, 165, 165];  
IVD.ys     = [136, 137, 135, 130, 130, 127, 123, 119, 117, 115, 112, 111, 109, 109, 111];      

IVD.N = length(IVD.zs);

% Z_L1_T12  = 22; %X_L1_T12  = 163; Y_L1_T12  = 112;

IVD.zs_interp = zeros(1,IVD.N);
for k = 1:IVD.N
    
z = IVD.zs(k);

% apply to slice
slice = get_slice_index(slices,levels,z,IVD.levels(k)) %#ok<NOPTS>

WI = WaterImages(:,:,slice);
binWI = BinaryWaterImages(:,:,slice);

% initial IVD mask seed
W = 11; % Width of intervertebral disc mask
H = 7; % Height of intervertebral disc mask 

% Middle of intevertabral disc coordinates
x = IVD.xs(k); 
y = IVD.ys(k); 

stop_flag = false;
stop_threshold = 50;
max_iteration  = 25;
iter = 1;
while ~stop_flag && (iter < max_iteration)
    [X, Y] = meshgrid(1:DW,1:DH); % Image size 
    EMask = ((X - x)/W).^2 + ((Y - y)/H).^2 <= 1;
    
BMask = (EMask.*binWI);
    mask_sum = sum(EMask(:)) - sum(BMask(:));
    if mask_sum < stop_threshold
        % expand mask
        W = W+1;
        H = H+1;
        % update the center of mass
        IVDmask = WI.*uint8(EMask).*uint8(binWI);
        CoM = center_of_mass(IVDmask);
        X = CoM(1);
        Y = CoM(2);
    else
        stop_flag = true;
    end
    iter = iter+1;
end

%% Apply to slice


figure,
imshow(WI)

IVDmask = WI.*uint8(EMask).*uint8(binWI);
figure,
imshow(IVDmask)

% Determine the intensities of the neighboring slices
stop_flag = false;
while ~stop_flag
   
   if slice > 2 && levels(slice-1) == levels(slice)
       WI_minus_1 = WaterImages(:,:,slice-1);
       binWI_minus_1 = BinaryWaterImages(:,:,slice-1);
       IVDmask_minus_1 = WI_minus_1.*uint8(EMask).*uint8(binWI_minus_1);
   else
       IVDmask_minus_1 = IVDmask;
   end
   if slice < (Nslices-1) && levels(slice+1) == levels(slice)
       WI_plus_1 = WaterImages(:,:,slice+1);
       binWI_plus_1 = BinaryWaterImages(:,:,slice+1);
       IVDmask_plus_1 = WI_plus_1.*uint8(EMask).*uint8(binWI_plus_1);       
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

end

IVD.zs_interp

%NOTE: add this to a IVD.zs array for all IVDs and 
% also construct a IVD.labels with this one labeled 12 for T12.
% The next one down will be labelled -1 for L1 and so forth.

if writeResults       
    save(resultsFile,'IVD');
end


