%-----------------------------------------------------------------------
% Author: Dr. Jason E. Hill, post-doctoral researcher
% Institution: Computational Neuroimaging Group at TTU
% Date: 1 JUL 2015
% Updated: 20 DEC 2016
%------------------------------------------------------------------------

clearvars, close all;
addpath(genpath(pwd))

% NOTE: For batch runs over many slices should set showProgress to 'false'
%       For exmining the analysis of an individual slice,
%       set showProgress to 'true' to see more infomation.
%       For individual slice run comment out the for i=1:Nslices line and
%       its 'end' line. Also comment out the saving of the .mat file.
%       Slices at the bottom of the heart cannot be viewed individually,
%       since they are masked from the slice above. Must run the for loop
%       from heartMax to the slice of interest.

%% choose subject and session

%% load data
%List of 4 intervention subjects: 1, 7, 8, 15
load subject07_i_registered.mat
resultsFile = 'subject07_i_results.mat';

SEdiamond3x3 = strel('diamond',1);
SEdiamond5x5 = strel('diamond',2);
SEdiamond7x7 = strel('diamond',3);
SEdiamond9x9 = strel('diamond',4);
SEdiamond11x11 = strel('diamond',5);
SEdiamond13x13 = strel('diamond',6);
SEdiamond15x15 = strel('diamond',7);
% image generation control
showProgress = false; %true;
batchShow    = ~showProgress;
writeResults = false;
oneBody = false; %true; % impose largest area for body mask (only do for subject 1)

niiFup_reslice = load_nii(NIFTI_file_name_F_upper_reslice);
niiFlo_reslice = load_nii(NIFTI_file_name_F_lower_reslice);
niiWup_reslice = load_nii(NIFTI_file_name_W_upper_reslice);
niiWlo_reslice = load_nii(NIFTI_file_name_W_lower_reslice);

%% slice loop inner

indeces = 1:Nslices;
% for article: subject 7, slices 44 and 116
% slice2show = 44
%indeces = 1:slice2show;
%indeces = bellybutton
%indeces = heartMax-5:heartApex+5
%indeces = (heartMax-1):(diaphram);%+15);
%indeces = (heartMax-1):(heartMax+1);
%indeces = 184:185
first_i  = indeces(1);
FH_i = length(upperSlices)+lowerTop-FH+1;

for i = indeces

if i == slice2show % for article: subject 15, slices 46 and 75
     showProgress = true;
     batchShow    = ~showProgress;
end
% problem slice:
%  if i == 61
%    showProgress = true;
%    batchShow    = ~showProgress;
%  end
i 

slice = slices(i);    
level = levels(i);    
voxelVolume = voxelVolumes(i);    

%% make 8-bit slice image of fat
% orient image slice (posterior down)
if level == 1
    IFraw = fliplr(rot90(niiFup_reslice.img(:,:,slice)));
elseif level == 2
    IFraw = fliplr(rot90(niiFlo_reslice.img(:,:,slice)));    
end
if showProgress && ~batchShow
    figure, imshow(IFraw, []);
end

IFrawSample = 0*IFraw;
IFrawSample(50:200,90:230) = IFraw(50:200,90:230);
IFrawSampleMax = max(max(IFrawSample));
IFraw(IFraw>IFrawSampleMax) = IFrawSampleMax;
% convert data to 8-bit gray-scale image
IFshow = uint8(double(IFraw)/double(IFrawSampleMax)*255);

if showProgress && ~batchShow
    figure, imshow(IFshow);
end

IFrawVeto = (IFshow < vetoFactorF*backgroundThreshold);
IFrawVeto = logical((1 - ~(IFshow < 0.75*vetoFactorF*backgroundThreshold).*imdilate(~IFrawVeto,ones(7))));
IFrawVeto = logical((1 - ~(IFshow < backgroundThreshold).*imdilate(~IFrawVeto,ones(9))));


%% make 8-bit slice image of water selected data

% orient image slice (posterior down)
if level == 1
    IWraw = fliplr(rot90(niiWup_reslice.img(:,:,slice)));    
elseif level == 2
    IWraw = fliplr(rot90(niiWlo_reslice.img(:,:,slice)));       
end

if showProgress && ~batchShow
    figure, imshow(IWraw,[]);
end
IWrawSample = 0*IWraw;
IWrawSample(50:200,90:230) = IWraw(50:200,90:230);
IWrawSample(110:135,135:175) = mean(mean(IWraw));
IWrawSampleMax = max(max(IWrawSample));
IWraw(IWraw > IWrawSampleMax) = IWrawSampleMax;
% convert data to 8-bit gray-scale image
IWshow = uint8(double(IWraw)/double(IWrawSampleMax)*255);

if showProgress && ~batchShow
    figure, imshow(IWshow);
end

%if level == 1
IWrawVeto = (IWshow < vetoFactorW*backgroundThreshold);
IWrawVeto = logical((1 - ~(IWshow < 0.75*vetoFactorW*backgroundThreshold).*imdilate(~IWrawVeto,ones(5))));
IWrawVeto = logical((1 - ~(IWshow < backgroundThreshold).*imdilate(~IWrawVeto,ones(3))));

    %else
%    IWrawVeto = (IWshow < 0.7*vetoFactorW*backgroundThreshold);
%end
% correct with IFveto mask
IallCorrect = imdilate(imfill(~IFrawVeto,'holes'),ones(3));
IWrawVeto = logical(~IallCorrect + IWrawVeto.*IallCorrect);

%% correct fat signal inhomogeneity

% construct rough estimate of the signal inhomogeneity via morphological closing
%if level == 1
[IFshow1,IFcorrect] = CorrectInhomogeneity(IFraw,34,[],1.0); %34 %100
%else
%    [IFshow1,IFcorrect] = CorrectInhomogeneity(IFraw,25); 
%end
IFshow1 = IFshow1.*uint8(1-IFrawVeto);
IFcorrect = IFcorrect.*(1-IFrawVeto);

if showProgress && ~batchShow
    figure, imshow(IFshow1);
end

%% segment the fat image foreground from background

% Find the peak and valley pixel bins of histogram
[peakPixelBinF,valleyPixelBinF] = FindHistogramExtrema(IFshow1,backgroundThreshold,foregroundThreshold1,(showProgress&&~batchShow));
%[peakPixelBinF,valleyPixelBinF] = FindHistogramExtrema(IFshow1.*uint8(1-IFrawVeto) + (IFshow.*uint8(IFrawVeto) + IFshow1.*uint8(IFrawVeto))/2,backgroundThreshold,foregroundThreshold1,(showProgress&&~batchShow));
peakPixelBinF
valleyPixelBinF

% rescale fat image using the histogram peak as the maximum  
peakIndicesF = find(IFshow1 == peakPixelBinF);
IFpeak = IFcorrect(peakIndicesF(1));
IFcorrect(IFcorrect>IFpeak) = IFpeak;
IFshow2 = uint8(double(IFcorrect)/double(IFpeak)*255);

if showProgress && ~batchShow
    figure, imshow(IFshow2);
end

% limit the valley pixel bin by Otsu's multithreshold method
valleyPixelBinF = AutorestrictThreshold(IFshow2,valleyPixelBinF,nbrThresholdsF);
valleyPixelBinF
% if valleyPixelBinF > 128
%     valleyPixelBinF  = valleyPixelBinF/2.0;
% end

% define fat foreground as being above histogram valley
seg_IFbin = (IFshow2 > valleyPixelBinF).*logical(1-IFrawVeto);

if showProgress && ~batchShow
    figure, imshow(seg_IFbin,[]);
end

%% correct water signal inhomogeneity

% construct rough estimate of the signal inhomogeneity via morphological closing
%if level == 1
    [IWshow2,~] = CorrectInhomogeneity(IWraw,50,[],1.0); %50 %100);
%else
%    [IWshow2,~] = CorrectInhomogeneity(IWraw,30,[],1.5);
%end

%% segment the image foreground from background within body cavity 

% mask corrected image with the internal (body cavity) mask
IWcorrect = double(IWshow2).*double(IWshow2>(IFshow2-4*backgroundThreshold));
%IWcorrect =
%double(IWshow2).*double(IWshow2>(IFshow2-6*backgroundThreshold)); with 15
%IWcorrect =
%double(IWshow2).*double(IWshow2>(IFshow2-4*backgroundThreshold)); with 10

% rescale image to the maximum intensity
IWmax = max(max(IWcorrect));
IWshow3 = uint8(double(IWcorrect)/double(IWmax)*255);
IWshow3 = IWshow3.*uint8(1-IWrawVeto);

if showProgress && ~batchShow
    figure, imshow(IWshow3);
end

% Find the peak and valley pixel bins of histogram
[peakPixelBinW,valleyPixelBinW] = FindHistogramExtrema(IWshow3,backgroundThreshold,foregroundThresholdW,(showProgress&&~batchShow));
%[peakPixelBinW,valleyPixelBinW] = FindHistogramExtrema(IWshow3.*uint8(1-IWrawVeto) + (IWshow.*uint8(IWrawVeto) + IWshow3.*uint8(IWrawVeto))/2,backgroundThreshold,foregroundThreshold1 - 2*backgroundThreshold,(showProgress&&~batchShow));
%[peakPixelBinW,valleyPixelBinW] = FindHistogramExtrema(IWshow3.*uint8(1-IWrawVeto) + (IWshow.*uint8(IWrawVeto) + IWshow3.*uint8(IWrawVeto))/2,backgroundThreshold,foregroundThreshold1 - 3*backgroundThreshold,(showProgress&&~batchShow));
%[peakPixelBinW,valleyPixelBinW] = FindHistogramExtrema(uint8((double(IWshow3).*((1+2.0*double(1-IWrawVeto))/3.0))),backgroundThreshold,foregroundThreshold1 - 2*backgroundThreshold,(showProgress&&~batchShow));
peakPixelBinW = peakPixelBinW - 2*backgroundThreshold
valleyPixelBinW = valleyPixelBinW %- backgroundThreshold

% rescale water image using the histogram peak as the maximum  
peakIndicesW = find(IWshow3 == peakPixelBinW);
if isempty(peakIndicesW)
    peakIndicesW = find(IWshow3 == peakPixelBinW + 1);
end

IWpeak = IWcorrect(peakIndicesW(1));
IWcorrect(IWcorrect>IWpeak) = IWpeak;
IWcorrect = IWcorrect.*(1-IWrawVeto);
IWshow4 = uint8(double(IWcorrect)/double(IWpeak)*255);

if showProgress && ~batchShow
    figure, imshow(IWshow4);
end

% save cross-section of saggital slice image
if i < Nslices
    Isaggital(:,Nslices-i+1) = IWshow3(:,160).*uint8(1-IWrawVeto(:,160));
end

% limit the valley pixel bin by Otsu's multithreshold method
valleyPixelBinW = AutorestrictThreshold(IWshow4,valleyPixelBinW,nbrThresholdsW);

% define water foreground as being above histogram valley
seg_IWbin = (IWshow4 > valleyPixelBinW).*double(1-IWrawVeto);

if showProgress && ~batchShow
    figure, imshow(seg_IWbin,[]);
end

% correct overlap with normalized raw signal
IallCorrect = logical((IFshow > IWshow).*(seg_IWbin.*seg_IFbin));
seg_IFbin(IallCorrect) = 1;
seg_IWbin(IallCorrect) = 0;
IallCorrect = logical(seg_IWbin.*seg_IFbin - IallCorrect);
seg_IFbin(IallCorrect) = 0;
seg_IWbin(IallCorrect) = 1;

%% separation of SCAT from VAT

% define a cut line mask to enable filling of holes apart from body cavity
cutMask = logical(0*seg_IFbin);
cutMask(size(seg_IFbin,1)-cutLineLength:size(seg_IFbin,1),round(0.5*size(seg_IFbin,2))) = true;
cutSaveLine = bwmorph((cutMask.*seg_IFbin),'clean');

% link SCAT together for belly button
linkMask = logical(0*seg_IFbin);
if sum(bellybutton == i) == 1
    linkMask(bbx,round(0.5*(size(seg_IFbin,2)-cutLineLength)):round(0.5*(size(seg_IFbin,2)+cutLineLength))) = true;
end  

% preliminary body and perimeter mask
bodyMask = imfill(logical(FindLargestArea(logical(1-IFrawVeto+linkMask+seg_IWbin))),'holes');
perimeterMask = logical(bwmorph(bodyMask,'dilate')-imerode(bodyMask,ones(7)));
% 
filledFatBase = logical(imfill((~cutMask.*(logical(seg_IFbin.*logical(1-seg_IWbin)+linkMask+perimeterMask))),'holes')+cutSaveLine);
% 
% % preliminary separation of SCAT from VAT using morphological openning 
openedBWfat = imdilate(imerode(filledFatBase,ones(5)),ones(5));  
% 
% % construct mask of entire body and perimeter more properly
bodyMask = imfill(bwmorph(bwmorph(logical(openedBWfat+linkMask+seg_IWbin+perimeterMask).*(seg_IWbin+seg_IFbin),'open'),'close'),'hole');
% construct mask of body perimeter
se1 = strel('rectangle',[11,11]);
perimeterMask = logical(bodyMask-imerode(bodyMask,se1));

% correct SCAT with inner cavity mask (and correct masks)
innerBase = imerode(imdilate(seg_IWbin.*(1-perimeterMask),ones(5)),ones(5)).*(1-linkMask);
innerBase(:,1:25) = 0;
innerBase(:,295:320) = 0;
% mask out the mammary tissue
if innerVeto(slice) > 0 && level == 1
    innerBase(1:innerVeto(slice),1:round(0.4*size(seg_IWbin,2)))      = 0;   %#ok<NASGU>
    innerBase(1:innerVeto(slice),round(0.6*size(seg_IWbin,2)):end)    = 0;   %#ok<NASGU>   
    innerBase(1:innerVeto(slice)+10,1:round(0.3*size(seg_IWbin,2)))   = 0;   %#ok<NASGU>
    innerBase(1:innerVeto(slice)+10,round(0.7*size(seg_IWbin,2)):end) = 0;   %#ok<NASGU>
    innerBase(1:innerVeto(slice)+20,1:round(0.25*size(seg_IWbin,2)))  = 0;   %#ok<NASGU>
    innerBase(1:innerVeto(slice)+20,round(0.75*size(seg_IWbin,2)):end)= 0;   %#ok<NASGU>
    innerBase(1:innerVeto(slice)+25,1:round(0.2*size(seg_IWbin,2)))   = 0;   %#ok<NASGU>
    innerBase(1:innerVeto(slice)+25,round(0.8*size(seg_IWbin,2)):end) = 0;   %#ok<NASGU>
    innerBase(1:innerVeto(slice)+30,1:round(0.15*size(seg_IWbin,2)))  = 0;   %#ok<NASGU>
    innerBase(1:innerVeto(slice)+30,round(0.85*size(seg_IWbin,2)):end)= 0;   %#ok<NASGU>  
end
innerAnterior = find(sum(innerBase,2)>0,1);
innerAnterior20 = find(sum(innerBase,2)>20,1);
if innerAnterior20 > (innerAnterior + 1)
    innerBase(1:innerAnterior20-1,:) = 0;
end
if strcmp(resultsFile,'subject08_f_results.mat') && (i > 95) && (i < 105)
    x = [1 round(size(innerMask,2)/2.0) 1 1];
    y = [round(size(innerMask,1)/2.0) size(innerMask,1) size(innerMask,1) round(size(innerMask,1)/2.0)];    
    IallCorrect  = logical(1-poly2mask(x,y,size(innerMask,1),size(innerMask,2)));
    innerBase = innerBase.*IallCorrect;
    x = [size(innerMask,2) round(size(innerMask,2)/2.0) size(innerMask,2) size(innerMask,2)];
    y = [round(size(innerMask,1)/2.0) size(innerMask,1) size(innerMask,1) round(size(innerMask,1)/2.0)];    
    IallCorrect  = logical(1-poly2mask(x,y,size(innerMask,1),size(innerMask,2)));
    innerBase = innerBase.*IallCorrect;
end
IallCorrect = fliplr(innerBase);
IallCorrect (:,1:round(size(innerBase,2)/3)) = 0;
IallCorrect (:,round(2*size(innerBase,2)/3):end) = 0;
% IallCorrect = IallCorrect + imerode(fliplr(innerBase),SEdiamond5x5);
% IallCorrect (:,1:round(0.2*size(innerBase,2))) = 0;
% IallCorrect (:,round(0.8*size(innerBase,2)):end) = 0;
innerBase = logical(innerBase + IallCorrect);
if showProgress && ~batchShow
    figure, imshow(innerBase,[]);
end
% frame the innerMask by the fat extrema
fatAnterior = find(sum(seg_IFbin,2)>0,1);
fatPosterior = size(seg_IFbin,1) - find(sum(flipud(seg_IFbin),2)>0,1) + 1;
fatLeftmost = find(sum(seg_IFbin,1)>0,1);
fatRightmost = size(seg_IFbin,2) - find(sum(flipud(seg_IFbin),1)>0,1) + 1;

innerBase(1:fatAnterior+10,:) = 0;
innerBase(fatPosterior-10:end,:) = 0;
innerBase(:,1:fatLeftmost+10) = 0;
innerBase(:,fatRightmost-10:end) = 0;

% further corrections and constructing of the inner cavity mask
innerBase = innerBase.*(1-linkMask);
innerBase = FindLargestArea(innerBase,'multiple',0.06*voidThreshold);
%innerBase = imerode(imfill(imdilate(innerBase,ones(9)),'holes'),SEdiamond9x9);
innerBase = imerode(imfill(imdilate(innerBase,ones(15)),'holes'),SEdiamond15x15);
innerBase = innerBase.*(1-linkMask);
IallCorrect = innerBase.*logical(logical((1-seg_IFbin)+seg_IWbin));
innerBase = logical(imerode(innerBase,ones(11))+IallCorrect);
innerBase = bwmorph(imfill(innerBase,'holes'),'close');
innerBase = FindLargestArea(innerBase,'multiple',0.1*voidThreshold);
if showProgress && ~batchShow
    figure, imshow(innerBase,[]);
end

%innerMask = logical(imdilate(imerode(innerBase,ones(5)),ones(5)).*logical(1-imdilate(perimeterMask,se1)));
innerMask = logical(imopen(innerBase,SEdiamond5x5).*logical(1-imdilate(perimeterMask,se1)));

% define the inner perimeter mask of the SCAT from the inner cavity mask
SCATinnerMask = bwmorph((imdilate(innerMask,se1)-innerMask),'close');
SCATinnerMaskFill = logical(bwmorph(imerode(imfill(imdilate(SCATinnerMask,ones(5)),'holes'),ones(5)),'close').*logical(1-imdilate(perimeterMask,se1)));
%SCATinnerMask = logical((SCATinnerMaskFill - imerode(SCATinnerMaskFill,se1)));
% correct inner mask
innerMask = logical(innerMask + imerode(SCATinnerMaskFill,se1));
innerMask = logical(logical(imdilate(innerMask,ones(5)).*innerBase) + innerMask);
innerMask = FindLargestArea(innerMask,'multiple',voidThreshold/5.0);
if showProgress && ~batchShow
    figure, imshow(innerMask,[]);
end

% connect any remaining gaps in the organs (bridge abdomenal muscle)
%if (innerVeto(slice) == 0 || level == 2) && i < (length(upperSlices)+lowerTop-FH-10)
% innerx = find((sum(imopen(imfill(imclose(innerBase,ones(11)),'holes'),ones(11)),2)==0)==0,1)+9;
% innerMask(innerx,round(0.5*(size(seg_IFbin,2)-1.75*cutLineLength)):round(0.5*(size(seg_IFbin,2)+1.75*cutLineLength))) = true;
% innerMask = logical(logical(imdilate(innerMask,ones(7)).*innerBase) + innerMask);
% innerMask = bwmorph(imerode(bwmorph(imfill(imdilate(innerMask,ones(3)),'holes'),'open'),ones(3)),'spur',10);
% %figure, imshow(innerMask,[]);
% innerMask = FindLargestArea(innerMask);
% if showProgress && ~batchShow
%     figure, imshow(innerMask,[]);
% end
SCATinnerMaskFill = innerMask;
if (i > (min(bellybutton)-2))&&(i < (max(bellybutton)-2))
    statsInner = regionprops(innerMask, 'ConvexHull');
    innerMask = roipoly(innerMask,statsInner(1).ConvexHull(:,1),statsInner(1).ConvexHull(:,2));
else
% correct inner mask with R & L convex hulls
% right side
IallCorrect = innerMask;

IallCorrect(1:round(size(innerMask,1)/2),1:round(0.42*size(innerMask,2)))   = 0;
if (i < FH_i - 50)
    IallCorrect(round(size(innerMask,1)/2)+1:end,:) = 0;
else
    IallCorrect(round(size(innerMask,1)/2)+1:end,1:round(size(innerMask,2)/2)) = 0;
end
NWpoint = find((IallCorrect(:,round(0.42*size(innerMask,2))+1)),1);
IallCorrect(NWpoint:round(size(innerMask,1)/2),round(0.42*size(innerMask,2))+1) = 1;
if showProgress && ~batchShow
    figure, imshow(IallCorrect,[]);
end
statsInner = regionprops(IallCorrect, 'ConvexHull');
SCATinnerMask = roipoly(IallCorrect,statsInner(1).ConvexHull(:,1),statsInner(1).ConvexHull(:,2));
if showProgress && ~batchShow
    figure, imshow(SCATinnerMask,[]);
end
% left side
IallCorrect = innerMask;
IallCorrect(1:round(size(innerMask,1)/2),round(0.58*size(innerMask,2)):end) = 0;
if (i < FH_i - 50)
   IallCorrect(round(size(innerMask,1)/2)+1:end,:) = 0;
else
   IallCorrect(round(size(innerMask,1)/2)+1:end,round(size(innerMask,2)/2)+1:end) = 0;
end

NWpoint = find((IallCorrect(:,round(0.58*size(innerMask,2))-1)),1);
IallCorrect(NWpoint:round(size(innerMask,1)/2),round(0.58*size(innerMask,2))-1) = 1;
IallCorrectLeftmost = find(IallCorrect(round(size(innerMask,1)/2),:)>0,1);
IallCorrect(round(size(innerMask,1)/2),IallCorrectLeftmost:round(0.58*size(innerMask,2))-1) = 1;
if showProgress && ~batchShow
    figure, imshow(IallCorrect,[]);
end
statsInner = regionprops(IallCorrect, 'ConvexHull');
IallCorrect = roipoly(IallCorrect,statsInner(1).ConvexHull(:,1),statsInner(1).ConvexHull(:,2));
if showProgress && ~batchShow
    figure, imshow(IallCorrect,[]);
end
innerMask = logical(IallCorrect + SCATinnerMask);
if showProgress && ~batchShow
    figure, imshow(innerMask,[]);
end
end

% correct inner mask with erosion and AT presence
IallCorrect = logical(seg_IFbin.*logical(innerMask - imerode(innerMask,SEdiamond3x3)));
innerMask = innerMask.*logical(1-IallCorrect);

% correct inner mask with spinal organ tissue    
%x = [round(size(innerMask,2)/3.0) round(size(innerMask,2)/2.0) round(2.0*size(innerMask,2)/3.0) round(size(innerMask,2)/3.0)];
%y = [size(innerMask,1) round(size(innerMask,1)/2.0) size(innerMask,1) size(innerMask,1)];
innerPosterior = size(SCATinnerMaskFill,1) - find(sum(flipud(SCATinnerMaskFill),2)>0,1) + 1;

% first correct the existing convex hull using the horzontal mid-line 
if (i < FH_i - 50)
IallCorrect = SCATinnerMaskFill;
IallCorrect(1:round(size(innerMask,1)/2.0),:) = false;
IallCorrect(innerPosterior-min(i/3,21)+1:end,:) = false;
IallCorrect((round(size(innerMask,1)/2.0))+1,:) = logical(IallCorrect((round(size(innerMask,1)/2.0))+2,:)+IallCorrect((round(size(innerMask,1)/2.0))+3,:));
IallCorrect = FindLargestArea(IallCorrect);
%IallCorrectLeftmost = find(sum((IallCorrect.*seg_IWbin),1)>0,1);
%IallCorrect(round(size(innerMask,1)/2.0)+1,IallCorrectLeftmost:round(size(innerMask,2)/2.0)) = true;
if showProgress && ~batchShow
    figure, imshow(IallCorrect,[]);
end
statsInner = regionprops(IallCorrect, 'ConvexHull');
IallCorrect = roipoly(IallCorrect,statsInner(1).ConvexHull(:,1),statsInner(1).ConvexHull(:,2));
if showProgress && ~batchShow
    figure, imshow(IallCorrect,[]);
end

% correct for any gross asymmetries
IallCorrect = logical(IallCorrect  + imerode(fliplr(IallCorrect),SEdiamond5x5));
statsInner = regionprops(FindLargestArea(SCATinnerMaskFill), 'ConvexHull');
SCATinnerMask = roipoly(SCATinnerMaskFill,statsInner(1).ConvexHull(:,1),statsInner(1).ConvexHull(:,2));
IallCorrect = IallCorrect.*SCATinnerMask; 
if showProgress && ~batchShow
    figure, imshow(IallCorrect,[]);
end
innerMask((round(size(innerMask,1)/2.0)+1):(innerPosterior-min(i/3,21)),:) = IallCorrect((round(size(innerMask,1)/2.0)+1):(innerPosterior-min(i/3,21)),:);
innerMask((round(size(innerMask,1)/2.0)),:) = logical(innerMask((round(size(innerMask,1)/2.0)),:)+innerMask((round(size(innerMask,1)/2.0)+1),:)+innerMask((round(size(innerMask,1)/2.0)-1),:));
innerMask(innerPosterior-min(i/3,21)+1:end,:) = SCATinnerMaskFill(innerPosterior-min(i/3,21)+1:end,:);
innerMask = imclose(innerMask,SEdiamond5x5);
if showProgress && ~batchShow
    figure, imshow(innerMask,[]);
end
end

% second built the spinal organ tissue mask
IallCorrect = 0*innerMask;
IallCorrect(innerPosterior-min(i/3,21):end,:) = true;

xS = [max(1,round(size(innerMask,2)/3.0)-i) round(size(innerMask,2)/2.0) min(size(innerMask,2),round(2.0*size(innerMask,2)/3.0)+i) max(1,round(size(innerMask,2)/3.0)-i)];
yS = [size(innerMask,1) round(size(innerMask,1)/2.0) size(innerMask,1) size(innerMask,1)];
IallCorrect = logical(IallCorrect + poly2mask(xS,yS,size(innerMask,1),size(innerMask,2)));

% correct inner mask with cleavage if needed 
if (i < heartMax-15) || (i > FH_i - 7) 
    xC = [1 round(size(innerMask,2)/2.0) size(innerMask,2) 1];
    yC = [1 round(size(innerMask,1)/2.0) 1 1];
    IallCorrect = logical(IallCorrect + poly2mask(xC,yC,size(innerMask,1),size(innerMask,2)));
end
IallCorrect = logical(IallCorrect.*FindLargestArea(bwmorph(imclose(seg_IWbin,SEdiamond5x5),'bridge'),'multiple',heartAreaThreshold)+(1-IallCorrect));

if showProgress && ~batchShow
    figure, imshow(IallCorrect,[]);
end
innerMask = imfill(innerMask.*IallCorrect,'holes'); 
if showProgress && ~batchShow
    figure, imshow(innerMask,[]);
end

% correct inner SCAT mask
SCATinnerMask = logical(imdilate(innerMask,se1)-innerMask);
%end
if ((i > diaphram + 1) || (i < FH_i - 50)) && ((i < (min(bellybutton)-2))||(i < (max(bellybutton)-2)))
innerAnterior = find(sum(innerMask,2)>0,1);
for k = 1:3
IallCorrect = SCATinnerMask;
IallCorrect(innerAnterior+16:end,:) = 0;
IallCorrect = imdilate(IallCorrect,SEdiamond3x3).*seg_IFbin;
SCATinnerMask = imdilate(SCATinnerMask,SEdiamond3x3).*seg_IFbin.*logical(1-innerMask);
if (sum(sum(IallCorrect(1:innerAnterior+15,:))) -  sum(sum(SCATinnerMask(1:innerAnterior+15,:)))) > 30;
    SCATinnerMask = logical(SCATinnerMask + IallCorrect);
    innerMask = innerMask.*logical(1-SCATinnerMask);
    SCATinnerMask = logical(imdilate(innerMask,se1)-innerMask);
end
end
end
SCATinnerMask = logical(imdilate(innerMask,se1)-innerMask);

openedBWfat = FindLargestArea(openedBWfat.*logical(1-innerMask)); % why problematic sometimes?

% apply cut line to perimeter SCAT ring
filledBWfat = logical(imfill((~cutMask.*logical(openedBWfat+perimeterMask+SCATinnerMask)),'holes')+cutSaveLine);

% accomplish more thorough separation of SCAT from VAT using erosion preserving surface
erodedBWfat = logical(imerode(filledBWfat,se1)+perimeterMask);

% clean up SCAT base mask
SCATbase = bwmorph(erodedBWfat,'open');
SCATbase = FindLargestArea(SCATbase);
if showProgress && ~batchShow
    figure, imshow(SCATbase,[]);
end

%% segment SCAT (subcutaneaous adipose tissue)

% remove any small spurs and indentations from SCAT base mask
SCATclosed = bwmorph(bwmorph(SCATbase,'spur'),'close');

% dilate base mask to encompass "true" area 
SCATdilate = imdilate(SCATbase,seSCAT);

% construct SCAT mask as the union of fat foreground with dilated base
% and removing any small spurs or isolated pixels
SCATmask = bwmorph(bwmorph(seg_IFbin.*SCATdilate,'clean'),'spur');
if oneBody
    SCATmask = FindLargestArea(SCATmask);
end

% define resulting SCAT image
SCAT = uint8(SCATmask).*IFshow2;

if ~batchShow
    figure, imshow(SCAT,[]);
end

% remove AT from water signal (organs)
seg_IWbin = seg_IWbin.*(1-seg_IFbin);

%% segment water bearing organs

% clean and fill holes in water foreground for organ mask
organsMask = bwmorph(bwmorph(logical(seg_IWbin),'clean'),'fill');

%% correct the SCAT mask

% define skin mask
seSkin = strel('rectangle',[5,7]);
skinMask = imdilate(logical(bodyMask-bwmorph(bodyMask,'erode')),se1);

% define cut line for SCAT
cutSaveLine = bwmorph((cutMask.*seg_IFbin),'clean');

% use inner boundary of body cavity mask to fill in SCAT properly
SCATfill = logical(imfill((logical(1-cutMask).*logical(SCATmask+SCATinnerMask)),'holes')+cutSaveLine).*seg_IFbin;
SCATcorrect = SCATfill.*logical(1-SCATmask);
SCATmask = logical(SCATmask + SCATcorrect);
 if ~batchShow
     figure, imshow(SCATmask,[]);
 end
SCATmask = FindLargestArea(SCATmask,'multiple',voidThreshold/50.0);
% if ~batchShow
%     figure, imshow(SCATmask,[]);
% end
SCAT = uint8(SCATmask).*IFshow2;

if ~batchShow
    figure, imshow(SCAT,[]);
end

% % correct segmentation of VAT 
% VATmask = seg_IFbin.*(1-SCATmask);

% % correct organs mask by removing voxels segmented to SCAT
% organsMask = organsMask.*logical(1 - SCATcorrect);

%% segment internal mask (body cavity)
% fill in the SCAT mask for an isolated thoracic body mask
se3 = strel('square',3);
se7 = strel('square',7);
SCATbase2 = logical(imerode(imclose(SCATmask,se7),se3)+perimeterMask);
SCATfill = imfill(SCATbase2,'holes');
if oneBody
    SCATfill = FindLargestArea(SCATfill);
end

% define the fill base mask for further void and VAT segmentation
fillbase = SCATfill-SCATbase2;
fillbase = FindLargestArea(fillbase);
if showProgress && ~batchShow
    figure, imshow(fillbase)
end
innerMask = imfill(logical(innerMask + imerode(fillbase,ones(5))),'holes');

% correct innerMask
SCATinnerMask = seg_IFbin.*logical(imdilate(innerMask,SEdiamond5x5)-innerMask);
innerMask = imclose(innerMask.*(1-SCATinnerMask),SEdiamond5x5);

% define the area outside the body cavity
externalMask = imopen(logical(1-innerMask),SEdiamond5x5);
externalMask(1,:) = 1;
externalMask(end,:) = 1;
externalMask(:,1) = 1;
externalMask(:,end) = 1;
if showProgress && ~batchShow
    figure, imshow(externalMask)
end
innerMask = logical(1-externalMask);
innerMask = imfill(innerMask,'holes');
externalMask = logical(1-innerMask);

% correct inner SCAT mask
SCATinnerMask = logical(imdilate(innerMask,se1)-innerMask);
SCATmask = SCATmask.*logical(1-innerMask);

SCATplusExternalMask = bwmorph(bwmorph(SCATmask,'clean'),'spur')|externalMask;
%SCATplusExternalMask = bwmorph(bwmorph(seg_IFbin.*SCATdilate,'clean'),'spur')|externalMask;
%figure, imshow(logical(1-SCATplusExternalMask),[]);
SCATmask = SCATmask.*SCATplusExternalMask;
% if ~batchShow
%     figure, imshow(SCATmask,[]);
% end
SCAT = uint8(SCATmask).*IFshow2;

%% segment voids (regions with low fat and water signal within the body cavity)

% define voids base as the regions not claimed by any existing mask
voidsBase = logical(1-seg_IWbin).*logical(1-organsMask).*logical(1-seg_IFbin).*logical(1-SCATplusExternalMask);

% define voids mask by opening its base
voidsMask = bwmorph(voidsBase,'open').*logical(1-organsMask).*logical(1-seg_IFbin).*logical(1-SCATplusExternalMask);
voidsMask(IFshow2 > min(valleyPixelBinW,valleyPixelBinF)) = 0;
voidsMask(IWshow4 > min(valleyPixelBinW,valleyPixelBinF)) = 0;
voids = uint8((peakPixelBinW/255.0)*(255-IWshow4).*uint8(voidsMask));

%% segment VAT (visceral adipose tissue)
% define VAT mask
VATmask = seg_IFbin.*(1-SCATmask);

if ~batchShow
    figure, imshow(VATmask,[]);
end

VATmask = VATmask.*innerMask.*logical(1-seg_IWbin);
VAT = uint8(VATmask).*IFshow2;
if ~batchShow
    figure, imshow(VAT,[]);
end

% define mask for non-fat areas
nonATmask = logical((1-SCATplusExternalMask).*(1-VATmask));
nonAT = uint8(nonATmask).*IFshow2;
if showProgress && ~batchShow
    figure, imshow(nonAT);
end


%% associate unassigned pixels

% needed for the remaining pixels not assigned to any mask
voidsCorrect = logical(voidsBase).*logical(1 - voidsMask);
% intensities of voxels near voids (fat saturated signal) 
NvoidsF = IFshow2.*uint8(voidsCorrect);
% intensities of voxels near voids (water saturated signal) 
NvoidsW = IWshow4.*uint8(voidsCorrect); %IWshow.*uint8(voidsCorrect);
% unassigned pixels are assigned based on the higher signal
organsCorrect = ((NvoidsW) > (NvoidsF)); 
% add unassigne voxels to appropriate masks
VATcorrect = logical((voidsCorrect - organsCorrect).*innerMask);
SCATcorrect = logical((voidsCorrect - organsCorrect - VATcorrect));
% update SCAT, VAT and organ masks
SCATmask(SCATcorrect) = 1;
SCAT = uint8(SCATmask).*IFshow2;
%figure, imshow(VATmask,[]);
VATmask(VATcorrect) = 1;
VAT = uint8(VATmask).*IFshow2;
organsMask(organsCorrect) = 1;
organs = IWshow.*uint8(organsMask);

if showProgress && ~batchShow
    figure, imshow(organs ,[]);
end

% Combine equally weighted data sets
%Iplus = sqrt((double(IFcorrect)/double(IFpeak)).^2 + (double(IWcorrect)/double(IWpeak)).^2);
Iplus = (IFcorrect/double(IFpeak)) + (IWcorrect/double(IWpeak));
IplusMax = max(max(Iplus));

% define "void" dominated image with in body cavity
IVraw = (IplusMax - Iplus);
if showProgress && ~batchShow
    figure, imshow(IVraw,[]);
end
IV = IVraw.*double(innerMask);
IV(isnan(IV)) = 0;
IV = (IV-1);
IV(IV<0) = 0;

bodyMask = logical(imfill(logical(1-(IVraw>1.5)),'holes') + imerode(bodyMask,ones(7)));
if oneBody && level == 1
    bodyMask(:,1:10) = 0;
    bodyMask(:,300:end) = 0;
    bodyMask = FindLargestArea(bodyMask);
end
SCATmask = bodyMask.*SCATmask;
organsMask = bodyMask.*organsMask;

%% segment mixed fractions

% create fat fraction image
FatFraction(:,:,i) = double(IFraw).*bodyMask./double(IFraw + IWraw);
if showProgress && ~batchShow
         figure, imshow(uint8(255*FatFraction(:,:,i)));
         colormap hot
end
FatFraction(:,:,i) = FatFraction(:,:,i).*(logical(1-IFrawVeto)|logical(1-IWrawVeto));

if showProgress && ~batchShow
         figure, imshow(FatFraction(:,:,i).*organsMask,[]);
end

% create water fraction image
WaterFraction(:,:,slice) = double(IWraw).*bodyMask./double(IFraw + IWraw);
if showProgress && ~batchShow
         figure, imshow(uint8(255*WaterFraction(:,:,slice)));
         colormap bone
end
WaterFraction(:,:,slice) = WaterFraction(:,:,slice).*(logical(1-IFrawVeto)|logical(1-IWrawVeto));
if showProgress && ~batchShow
         figure, imshow(WaterFraction(:,:,slice).*VATmask,[]);
end


%% aorta segmentation

% define the aorta base mask from seed
aortaBase = 0.*organsMask;
if level == 1
%    if (i > first_i) && (i > 15) && (sum(sum(aortaSeed(:,:,slice))) == 0) && (sum(sum(PAAT3d(:,:,i-1)))>0)
%        aortaSeed(:,:,slice) = aortaSeed(:,:,slice+1);
%    end   
    aortaBase = aortaSeed(:,:,slice);
    % take care of continuity

end
    
% grow aortaMask from seed by dilation
sea = strel('disk',aortaRadius-round(i/25));
aortaMask = imdilate(aortaBase,sea).*(IWshow>aortaThreshold);
aortaMask = imfill(aortaMask,'holes');

% choose largest areas within field as aorta near the image center
[labeledBWaorta,nbr_labels_aorta] = bwlabel(aortaMask,4);
statsAorta = regionprops(labeledBWaorta,'area','centroid');
label_aorta_areas = cat(1, statsAorta.Area);
label_aorta_centroids = cat(1, statsAorta.Centroid);
k_aorta_max = max(label_aorta_areas);
aorta_maxLabel = find(label_aorta_areas == k_aorta_max);
% selection of largest area is restricted by certain conditions
if ~isempty(k_aorta_max)
    if k_aorta_max > aortaMinArea && ((label_aorta_centroids(aorta_maxLabel,1)-round(0.5*size(aortaMask,1))...
            + label_aorta_centroids(aorta_maxLabel,2)-round(0.5*size(aortaMask,2)))<aortaDistance);
        aortaMask = bwmorph((labeledBWaorta == aorta_maxLabel),'open');
    else
        aortaMask = 0*aortaMask;
    end
else
    aortaMask = 0*aortaMask;
    PAATmask = aortaMask;
    PAAT = IFshow2.*uint8(PAATmask);    
end
label_aorta_areas2 = find(label_aorta_areas > 0.5*k_aorta_max & label_aorta_areas < k_aorta_max);
if ~isempty(label_aorta_areas2)
    aorta_maxLabel = label_aorta_areas2(1);
    k_aorta_max = label_aorta_areas(aorta_maxLabel);
    if ~isempty(k_aorta_max)
        if k_aorta_max > aortaMinArea && ((label_aorta_centroids(aorta_maxLabel,1)-round(0.5*size(aortaMask,1))...
             + label_aorta_centroids(aorta_maxLabel,2)-round(0.5*size(aortaMask,2)))<aortaDistance);
         aortaMask = logical(aortaMask + bwmorph((labeledBWaorta == aorta_maxLabel),'open'));
        end
    end
end

aortaMask = logical(aortaMask + aortaBase); 
aorta = IWshow4.*uint8(aortaMask);

% define aorta veto line for heart segmentation from aorta seed
aortaVeto = aortaBase;
sum(sum(aortaBase))
if sum(sum(aortaBase)) < 3
    aortaMax = find(sum(aortaMask,2)>0,1);
    aortaVeto(aortaMax+6:end,:) = 1;                
%    aortaVeto = imdilate(aortaBase,ones(5,50));
else
    % curve of aorta, no heart
    aortaBase = aortaSeed(:,:,slice-3); 
    aortaVeto = logical(1 + aortaMask);
end
aortaAnterior = find(sum(aortaMask,2)>0,1);
%aortaAnterior = find(sum(aortaBase,2)>0,1);
if isempty(aortaAnterior)
    aortaAnterior = round(size(innerMask,1)/2.0);
end


%% Lung segmentation
if i < diaphram + 21    % allow lungs to extend a bit past the diaphram

%% segment PAAT (periaortic adipose tissue)

% define the PAAT base by using the aortaMask as a seed
PAATbase = aortaMask;
for n = 1:PAATdilate
    PAATbase = logical((bwmorph(imdilate(PAATbase,SEdiamond3x3),'bridge').*VATmask)|aortaMask);
end
if showProgress && ~batchShow
    figure, imshow(PAATbase,[]);
end
% take care of continuity
if (i > first_i) && (i > 15)
    PAATbase = PAATbase + (PAAT3d(:,:,i-1)>0);
end

%PAATbase = bwmorph(bwmorph(aortaFatBase,'dilate'),'bridge');
if PAATbase > voidThreshold/5.0
    PAATmask = (FindLargestArea(PAATbase).*VATmask);
else
    PAATmask = PAATbase.*VATmask;
end
   

% correct PAATmask with aorta (require to not be above the aorta)
% sePAATcorrect = ones(1+2*PAATdilate);
% sePAATcorrect(1:PAATdilate,:) = 0;
% PAATcorrect = imdilate(aortaMask,sePAATcorrect);
PAATcorrect = PAATmask;
PAATcorrect(:) = true;
PAATcorrect(1:(aortaAnterior-1),:) = false;
PAATcorrect = logical(PAATcorrect + bwmorph(aortaMask,'dilate',3));
PAATmask = logical(PAATmask.*PAATcorrect);
PAAT = IFshow2.*uint8(PAATmask);

PAATanterior = find(sum(PAATmask,2)>0,1);
PAATposterior = size(PAATmask,1) - find(sum(flipud(PAATmask),2)>0,1) + 1;    
 
% correct the aorta mask
aortaMask = imclose(PAATmask + aortaMask,SEdiamond3x3).*organsMask;
aorta = IWshow4.*uint8(aortaMask);

%% correct void signal inhomogeneity

% construct rough estimate of the signal inhomogeneity via morphological closing
[IVshow,~] = CorrectInhomogeneity(IV,150,'medfilt',1.0);

if showProgress && ~batchShow
    figure, imshow(IVshow,[]);
end

% find peak and valley in histogram of the void image foreground
%foregroundThresholdV = 128;
%[peakPixelBinV,valleyPixelBinV] = FindHistogramExtrema(IVshow,backgroundThreshold,foregroundThresholdV,(showProgress&&~batchShow));
peakPixelBinV = 255
valleyPixelBinV = 128
%valleyPixelBinV = 128
% adjust valley back to original value
% valleyPixelBinV = round(valleyPixelBinV*peakPixelBinV/255);
% valleyPixelBinV

% restrict foreground threshold to be equal to or less than the next to
% highest automatically determined threshold
valleyPixelBinV = AutorestrictThreshold(IVshow,valleyPixelBinV,nbrThresholdsV);
valleyPixelBinV

% define void foreground base as being above histogram valley, use for lungs
lungBase = bwmorph((IVshow > valleyPixelBinV),'open');
lungBase = imclose(lungBase,ones(11)).*bwmorph((IVshow > valleyPixelBinV-4.0*backgroundThreshold),'open');
%lungBase = logical(lungBase.*(1-seg_IFbin).*(1-seg_IWbin));
lungBase = imfill(imclose(lungBase,ones(3)),'holes');
if showProgress && ~batchShow
    figure, imshow(lungBase,[]);
end

% lungMask = FindLargestArea(lungBase,'multiple',0.04*voidThreshold);
% %lungMask = FindLargestArea(lungBase,'multiple',0.1*voidThreshold);
% lungMask = imclose(lungMask,ones(17));
% IallCorrect = logical((IFshow.*uint8(lungMask) < 1.25*vetoFactorF*backgroundThreshold).*seg_IFbin.*lungMask);
% seg_IFbin(IallCorrect) = 0;
% lungBase(IallCorrect)  = true;
% IallCorrect = logical((IWshow.*uint8(lungMask) < 2.0*vetoFactorW*backgroundThreshold).*seg_IWbin.*lungMask);
% seg_IWbin(IallCorrect) = 0;
% lungBase(IallCorrect)  = true;

% handle points within the lungs
lungMask = logical(imerode(FindLargestArea(lungBase,'multiple',0.04*voidThreshold),ones(7))+(IVshow == 255));
lungMask = imclose(lungMask,ones(3));
if showProgress && ~batchShow
     figure, imshow(lungMask,[]);
end
lungBase = logical(lungBase.*(1-seg_IFbin).*(1-seg_IWbin));
IallCorrect = logical((IFshow.*uint8(lungMask) < 1.25*vetoFactorF*backgroundThreshold).*seg_IFbin.*lungMask);
seg_IFbin(IallCorrect) = 0;
lungBase(IallCorrect)  = true;
IallCorrect = logical((IWshow.*uint8(lungMask) < 2.5*vetoFactorW*backgroundThreshold).*seg_IWbin.*lungMask);
seg_IWbin(IallCorrect) = 0;
lungBase(IallCorrect)  = true;

% handle points on the edge of the lungs
lungMask = lungBase.*logical(1-lungMask);
if showProgress && ~batchShow
     figure, imshow(lungMask,[]);
end
IallCorrect = logical((IFshow2.*uint8(lungMask) > valleyPixelBinF - 2.0*backgroundThreshold));
seg_IFbin(IallCorrect) = true;
lungMask(IallCorrect) = 0;
IallCorrect = logical((IWshow4.*uint8(lungMask) > valleyPixelBinW - 2.0*backgroundThreshold));
seg_IWbin(IallCorrect) = true;
lungMask(IallCorrect) = 0;
lungBase = logical((lungBase + lungMask).*(1-seg_IFbin).*(1-seg_IWbin).*(1-PAATmask));
lungBase = logical(FindLargestArea(lungBase,'multiple',0.04*voidThreshold));
lungBase = imclose(lungBase,ones(3));
seg_IFbin(lungBase) = 0;
seg_IWbin(lungBase) = 0;

% define lungMask
if i < diaphram + 1
    % assume largest voids are lung tissue
    lungMask = logical(logical(FindLargestArea(lungBase,'multiple',0.04*voidThreshold)).*logical(1-PAATmask));
else
    % assume no lung tissue
    lungMask = 0.*lungBase;
end
if i > first_i+1     
%    figure, imshow(lungBase,[]);
    lungMask = logical(logical(lungMask + logical(lungBase+voidsMask).*imdilate((lung3d(:,:,i-1)>0),SEdiamond3x3)));
end
if showProgress && ~batchShow
     figure, imshow(lungMask,[]);
end

%lungMask = logical(FindLargestArea(lungBase,'multiple',0.1*voidThreshold));

% correct lung mask with nearby voids
for n = 1:3
   lungMask = logical(imdilate(lungMask,SEdiamond3x3).*logical(voidsMask+lungMask));
end
%lungMask = bwmorph(FindLargestArea(lungBase,'multiple',0.01*voidThreshold),'close');

% correct existing masks
voidsCorrect = logical(uint8(voidsMask-lungMask));
voidsMask(lungMask) = 0;
voidsMask(PAATmask) = 0;
%figure, imshow(VATmask,[]);
VATcorrect = ((VAT.*uint8(voidsCorrect))>=(organs.*uint8(voidsCorrect)))&voidsCorrect;
VATmask(VATcorrect) = 1;
VATmask = VATmask.*logical(1-lungMask).*logical(1-voidsMask);
PAATmask = PAATmask.*logical(1-lungMask).*logical(1-voidsMask);
%figure, imshow(VATmask,[]);
organsCorrect = logical(logical(voidsCorrect).*logical(1 - VATcorrect));
organsMask(organsCorrect) = 1;
organsMask = organsMask.*(1-lungMask);


%% heart segmentation

% find most anterior point of lung
%lungAnterior = find(sum(lungMask,2)>0,1);
lungAnterior = min([find(sum(lungMask,2)>0,1),find(sum(voidsMask,2)>0,1)]);

% define left lung veto mask to aid in heart segmentation
leftLungVetoMask = IVshow;
if (i > heartMax) && (i < (heartApex + 1))
    leftLungVetoMask = leftLungVetoMask + 100*uint8(bwmorph(imdilate(organsMask,ones(3)),'perim4'));
    leftLungVetoMask(:,1:round(size(IVshow,2)/2)+10) = 0;
    leftLungVetoMask = FindLargestArea(imopen(imfill(imclose((leftLungVetoMask > 90),ones(7)),'holes'),ones(3)));
    leftLungVetoMask = logical(leftLungVetoMask + imdilate(leftLungVetoMask,ones(7)).*(IVshow > 74));
    leftLungVetoMask = imerode(imfill(imdilate(leftLungVetoMask,ones(5)),'holes'),ones(5));
%    figure, imshow(leftLungVetoMask,[]);
    IallCorrect = logical(1 - imdilate(organsMask,ones(13)).*imdilate(IVshow > 90,ones(13)));
    IallCorrect(lungAnterior + 30:end,:) = true;
%    figure, imshow(IallCorrect,[]);
    IallCorrect = IallCorrect.*imfill(imdilate(IVshow > 74,ones(11)),'holes');
%    figure, imshow(IallCorrect,[]);
    leftLungVetoMask = leftLungVetoMask.*IallCorrect;
    leftLungVetoMask = FindLargestArea(leftLungVetoMask,'multiple',0.2*voidThreshold);
else
    leftLungVetoMask = lungMask;  
    leftLungVetoMask(:,1:round(size(IVshow,2)/2)+5) = 0;
end
if showProgress && ~batchShow
    figure, imshow(leftLungVetoMask,[]);
end
% % define right lung veto mask to aid in heart segmentation
% rightLungVetoMask = IVshow;
% if (i > heartMax) && (i < (heartApex + 1))
%     rightLungVetoMask = rightLungVetoMask + 100*uint8(bwmorph(imdilate(organsMask,ones(3)),'perim4'));
%     rightLungVetoMask(:,round(size(IVshow,2)/2)-5:end) = 0;
%     rightLungVetoMask = FindLargestArea(imopen(imfill(imclose((rightLungVetoMask > 90),ones(7)),'holes'),ones(3)));
%     IallCorrect(lungAnterior + 20:end,:) = true;
%     rightLungVetoMask = rightLungVetoMask.*IallCorrect;
%     rightLungVetoMask = FindLargestArea(rightLungVetoMask,'multiple',0.2*voidThreshold);
% else
%     rightLungVetoMask = lungMask; 
%     rightLungVetoMask(:,round(size(IVshow,2)/2)-5:end) = 0;
% end    

% restrict to only above observed bottom of heart
if i > heartApex
    % do nothing
    heartMask = 0*organsMask;
else
    % assume that heart is near the image center (shifted up and to the left side of body)
    heartShift = min((2*i),heartShift);   
    heartBase = 0*organsMask;
    heartBase(round(0.5*(size(organsMask,1)-heartShift)),round(0.5*size(organsMask,2))+heartShift) = 1;
    if showProgress && ~batchShow
        figure, imshow(heartBase,[]);
    end
    seh = strel('disk',heartDilate);
    %heartBase = imfill(heartBase,'holes');
    heartBase = imdilate(heartBase,seh).*logical(organsMask + VATmask).*logical(1-aortaVeto).*logical(1-leftLungVetoMask);
    heartBase = imfill(heartBase,'holes');
    heartBase(1:lungAnterior-1,:) = 0;
    if showProgress && ~batchShow
        figure, imshow(heartBase,[]);
    end

    % use lungs to frame the heart region
    seNotlung = strel('disk',lungErode);
    notLung = logical((1-lungMask).*fliplr(1-lungMask));
    %     figure, imshow(lungMask,[])
    %         figure, imshow(heartBase,[])
    notLung = imerode(notLung.*heartBase,seNotlung);
    if sum(sum(notLung)) < 400
        notLung = FindLargestArea(notLung);
    end

    % figure, imshow(notLung,[])

    % plot lung mask and centroid framed by lung
    statsHeart = regionprops(notLung,'centroid');
    centroids = cat(1, statsHeart.Centroid);

    if showProgress && ~batchShow && ~isempty(centroids)
        figure, imshow(lungMask)
        hold on
        plot(centroids(:,1),centroids(:,2), 'b*')
        hold off
    end

    centroidMask = 0*notLung;
    if ~isempty(centroids)
        centroidMask(round(centroids(:,2)),round(centroids(:,1))) = 1;
    end

    % remove side lobes that can be caused by diaphraim ...
    lungVeto = bwmorph(imfill(imclose(lungMask,ones(35)),'holes')...
                   -imclose(lungMask,ones(35)),'dilate',6)...
                   .*logical(1-lungMask);
    % ... unless it is the heart region.              
    if sum(sum(logical(lungVeto - lungVeto.*logical(1-centroidMask)))) == 1
        lungVeto = 0.*lungVeto;
    end

    % assume that the heart is framed by the lungs
    notLung = FindLargestArea(notLung);
    notLung = bwmorph(notLung.*logical(1-lungVeto),'open');

    if showProgress && ~batchShow
        figure, imshow(notLung,[]);
    end
    %figure, imshow(logical(1-lungMask),[])
    %figure, imshow(logical(organsMask+VATmask),[])

    % update heart base by dilation and masking
    seHeart = strel('disk',notLungDilate);  %21
    if (i > heartMax) && (i < (heartApex + 1))
        seInner = strel('disk',round(0.5*innerDilate));    %13    
    else
        seInner = strel('disk',innerDilate);    %13
    end
    %seInnerVeto = strel('disk',round(0.5*innerDilate));    %13
    %SCATfill = logical(logical(1-imdilate(SCATinnerMask,seInner))+imdilate(lungMask,seInner));
    heartBase = imdilate(notLung,seHeart).*logical(1-lungMask)...
               .*logical(organsMask+VATmask)...
               .*logical(1-bwmorph(aortaMask,'dilate'))...
               .*logical(1-imdilate(SCATinnerMask,seInner))...
               .*logical(1-aortaVeto);
    % get rid of spurs
    %figure, imshow(heartBase,[]);
    heartBase = imclose(heartBase,ones(5));
    %heartBase = imerode(heartBase,ones(5));
    %figure, imshow(heartBase,[]);
    if showProgress && ~batchShow
        figure, imshow(heartBase,[]);
    end

    % define the heart mask
    heartMask = imdilate(FindLargestArea(heartBase),ones(5));
    %figure, imshow(heartMask,[]);
    heartMask = bwmorph(heartMask.*logical(1-lungVeto),'open').*organsMask.*logical(1-aortaVeto);
    if i < (heartApex - 1)
        heartMask = FindLargestArea(heartMask,'multiple',heartAreaThreshold);
    else
        heartMask = FindLargestArea(heartMask,'multiple',round(heartAreaThreshold/10));
    end
    if showProgress && ~batchShow
        figure, imshow(heartMask,[]);
    end

    % mask lower part of heart with previous segmentation to restrict its area
    if i == heartMax
        heartMaxLeft = find(sum(heartMask)>0,1);
    end

    if (i > heartMax) && (i < (heartApex + 1)) && (i > first_i+1)
    %      if showProgress && ~batchShow
    %          figure, imshow(heartMask,[]);
    %      end
    %    if i < (heartApex-2)
         seHC = strel('disk',1);    
    %      seHC = strel('disk',CATdilate-3);
    %    else
    %      seHC = strel('disk',CATdilate-4);
    %    end
        heartCorrect = logical(heartMask.*imdilate(logical(imerode((heart3d(:,:,i-1)>0)+(EAT3d(:,:,i-1)>0),ones(5))),ones(7)).*logical(1-leftLungVetoMask)); 
     %    heartCorrect = logical(heartMask.*imdilate(logical(imerode((heart3d(:,:,i-1)>0)+(EAT3d(:,:,i-1)>0),ones(5 + floor(0.5*(i - heartMax))))),ones(7 + floor(0.5*(i - heartMax)))).*logical(1-leftLungVetoMask)); 
        % mask heart with the aorta fat, moving up as approach apex
        heartCorrect(PAATanterior - 3.0*(i - heartMax):end,:) = 0;
        heartCorrect(:,1:heartMaxLeft+(i - heartMax)) = 0;
        heartMean = mean(double(IWshow(logical(heartCorrect))));
        heartStd = std(double(IWshow(logical(heartCorrect))));
         if showProgress && ~batchShow
             figure, imshow(heartCorrect,[]);
         end
    
        diaphramCorrect = heartMask.*logical(1-heartCorrect).*logical(1-leftLungVetoMask).*logical(1-imdilate(imdilate(SCATinnerMask,seInner),SEdiamond7x7));
        diaphramCorrect = FindLargestArea(diaphramCorrect);
        if showProgress && ~batchShow
             figure, imshow(diaphramCorrect,[]);
        end
        IallCorrect = logical(imdilate(logical(diaphramCorrect),ones(17)).*organsMask.*(1-heartCorrect));
        diaphramSum = sum(sum(diaphramCorrect));
        diaphramMean = mean(double(IWshow(logical(diaphramCorrect))));
        diaphramStd = std(double(IWshow(logical(IallCorrect))));    
        diaphramStd = max([heartStd diaphramStd]);
        if diaphramSum > 64 && diaphramMean > (heartMean + 5.0)
            heartCorrectThresh = heartMean + 0.5*(diaphramMean-heartMean)*(heartStd/diaphramStd);
    %        heartCorrectThresh = 0.5*(heartMean + heartStd + diaphramMean - diaphramStd);
    %        heartCorrectThresh = 0.5*(heartMean + diaphramMean);
            if  heartCorrectThresh > heartMean
                heartCorrect = bwmorph(logical(heartCorrect).*((IWshow.*uint8(heartCorrect))<heartCorrectThresh),'close');
                heartCorrectFlag = true;  
            else
                heartCorrectFlag = false; 
            end
        end   
    %    heartCorrect = logical(heartCorrect + heartMask.*imdilate(VATmask,ones(7)) + heartMask.*imdilate(lungMask,ones(5)));
    %    heartCorrect = imdilate(imerode(heartCorrect,ones(5)),ones(5));
        heartCorrect = imerode(imdilate(heartCorrect,ones(6)),ones(5));
        heartCorrect = logical(imfill(heartCorrect+imdilate(heartCorrect,ones(9)).*VATmask,'holes'));
        heartMask = organsMask.*heartCorrect;
        heartMask  = FindLargestArea(heartMask);
        if showProgress && ~batchShow
            figure, imshow(heartMask,[]);
        end  
    % else
    %     heartCorrect = imdilate(imerode(heartCorrect,ones(5)),ones(5));
    %     heartMask = organsMask.*heartCorrect;
    end
    
    %figure, imshow(heartMask,[])

    % mask heart from lung extreme
    heartMask(1:lungAnterior-1,:) = 0;
    % mask heart from tissue near lungs/voids
    heartMask = FindLargestArea(organsMask.*imdilate(heartMask.*(1-imdilate(logical(lungMask+voidsMask+aortaMask+PAATmask+VATmask+imdilate(SCATinnerMask,seInner)),SEdiamond7x7)),SEdiamond7x7),'multiple',round(heartAreaThreshold/10));
end

if showProgress && ~batchShow
   figure, imshow(heartMask,[]);
end  
sum(heartMask(:))
%figure, imshow(heartMask,[])

%% segment CAT (any thoracic adipose tissue near heart) AKA Peri-cardial Adipose Tissue 

% define CAT via dilation
if (i > heartMax) && (i < (heartApex + 1)) && (i > first_i+1)
%    seCAT = strel('disk',CATdilate);
%    seInner2 = strel('disk',innerDilate-CATdilate);
%    heartFatBase = bwmorph(bwmorph((imdilate(heartCorrect,seCAT).*logical(1-imdilate(SCATinnerMask,seInner2)).*VATmask)|heartCorrect,'spur'),'open');;
    heartFatBase = logical(heartCorrect + imdilate(FindLargestArea(heartBase),ones(5)));
else
    seCAT = strel('disk',CATdilate);% + CATdilateCorrect(i));
%    innerDilate2 = max([innerDilate-CATdilate,CATdilate-innerDilate,1]);
%    seInner2 = strel('disk',innerDilate2);
%    heartFatBase = bwmorph(bwmorph((imdilate(heartMask,seCAT).*logical(1-imdilate(SCATinnerMask,seInner2)).*VATmask)|heartMask,'spur'),'open');
    heartFatBase = bwmorph(bwmorph((imdilate(heartMask,seCAT).*logical(1-imdilate(SCATinnerMask,ones(3))).*VATmask)|heartMask,'spur'),'open');
%    CATbase = bwmorph(bwmorph(heartFatBase,'dilate'),'bridge');
end
CATbase = imfill(bwmorph(imdilate(heartFatBase,SEdiamond5x5),'bridge'),'holes').*logical(1-leftLungVetoMask); %.*logical(1-rightLungVetoMask);
CATbase(1:lungAnterior-7,:) = 0;  

% find largest fat areas
CATmask = (FindLargestArea(logical(CATbase.*VATmask + heartMask)).*VATmask); 
if i > (heartMax)
    CATmask = FindLargestArea(CATmask,'multiple',CATareaThreshold);
else
    CATmask = FindLargestArea(CATmask,'multiple',round(CATareaThreshold/10.0));
end

if showProgress && ~batchShow
   figure, imshow(CATmask,[]);
end

%if i < (heartMax+1) || ~heartCorrectFlag
% correct heart mask with CAT mask
CATleftmost = find(sum(CATmask)>0,1);
%CATmaskSum = logical(sum(CATmask));
%CATlen = length(CATmaskSum);
%CATrightmost = CATlen;
% for j = CATlen:-1:1
%     if CATmaskSum(j) && (CATrightmost == CATlen)
%         CATrightmost = j;
%     end
% end
CATleftmost = find(sum(CATmask)>0,1);
if ~isempty(CATleftmost)
CATveto = CATmask;
if CATleftmost > CATmargin + 15
    CATveto(:,CATleftmost-CATmargin - 15:size(CATbase,2)) = 1;
else
    CATveto = ones(size(CATmask,1),size(CATmask,2));
end
heartMask = heartMask.*CATveto;
heartMask = FindLargestArea(heartMask,'multiple',heartAreaThreshold);
end
%CATmask = logical(CATmask + imdilate(heartMask,ones(5)).*VATmask);
%heartMask = imclose(logical(CATmask+heartMask),ones(5)).*organsMask;

heartPosterior = size(heartMask,1) - find(sum(flipud(heartMask),2)>0,1) + 1;
CATmask(heartPosterior+CATmargin:size(CATbase,1),:) = 0;

%end

% % define edge mask of void image to attempt to restrict to pericardial sack
% CATsat = IFshow2.*uint8(bwmorph(CATmask,'erode'));
% CATsat(CATsat==0) = 255;
% CATedgeBase = bwmorph(128-CATsat,'clean');
% CATedgeBase = CATedgeBase.*(1-imdilate(lungMask,ones(9)));
% CATedgeBase = bwmorph(CATedgeBase,'clean');
% seCEM = [[1,1,1,1,1];[0,1, 1, 1,0];[0,0,1,0,0];[0,0,0,0,0];[0,0,0,0,0]];
% CATedgeMask = imdilate(CATedgeBase,seCEM);
% CATmask = CATmask.*(1-CATedgeMask);

% fill in any holes in CAT %and heart mask
CATheart = imfill(bwmorph(bwmorph(logical(CATmask+heartMask),'dilate',2),'erode',2),'holes');
%heartMask = logical(heartMask + organsMask.*CATheart);
CATmask = logical(CATmask + VATmask.*CATheart);

% correct CAT mask
CATmask = CATmask.*logical(1-PAATmask);


%% associate unassigned pixels

% needed for the remaining pixels not assigned to any mask
voidsCorrect = logical(bwmorph(bodyMask,'erode',3)).*logical(1 - (SCATmask + VATmask + organsMask + voidsMask + CATmask + PAATmask + lungMask + heartMask + aortaMask));
% intensities of voxels near voids (fat saturated signal) 
NvoidsF = IFshow2.*uint8(voidsCorrect);
% intensities of voxels near voids (water saturated signal) 
NvoidsW = IWshow4.*uint8(voidsCorrect);
% intensities of voxels near voids (void saturated signal) 
NvoidsV = IVshow.*uint8(voidsCorrect);
% unassigned pixels are assigned based on the higher signal
organsCorrect = ((NvoidsW) > (NvoidsF))&((NvoidsW) > (NvoidsV)); 
fatCorrect = ((NvoidsF) > (NvoidsW))&((NvoidsF) > (NvoidsV)); 
voidsCorrect = logical(logical(voidsCorrect).*logical(1 - organsCorrect).*logical(1 - fatCorrect).*logical(1-PAATmask));
% add unassigne voxels to appropriate masks
VATcorrect = logical((fatCorrect).*innerMask);
SCATcorrect = logical(logical(fatCorrect).*logical(1 - VATcorrect));
% update SCAT, VAT and organ masks
SCATmask(SCATcorrect) = 1;
SCAT = uint8(SCATmask).*IFshow2;
VATmask(VATcorrect) = 1;
VAT = uint8(VATmask).*IFshow2;
organsMask(organsCorrect) = 1;
organs = IWshow4.*uint8(organsMask);
voidsMask(organsCorrect) = 0;
voidsMask(fatCorrect) = 0;
voidsMask(voidsCorrect) = 1;
voids = IVshow.*uint8(voidsMask);

% adipose tissue array
lungs    = IVshow.*uint8(lungMask);
ATarray  = [SCAT,VAT,organs,lungs];

tissues(:,:,1) = 0.5*organs + SCAT;
tissues(:,:,2) = VAT + SCAT;
tissues(:,:,3) = lungs + 0.5*voids;

if showProgress && ~batchShow
    figure, imshow(tissues);
end


%% segment EAT (epicardial adipose tissue)
% preprocess heartMask to exclude spurious tissue near lungs
statsEAT = regionprops(heartMask, 'ConvexHull');
if isempty(statsEAT)
    EATmask = 0.*heartMask;
else
    EATmask = roipoly(heartMask,statsEAT(1).ConvexHull(:,1),statsEAT(1).ConvexHull(:,2)).*VATmask;
    CATmask = logical(CATmask + EATmask);
end
% IallCorrect = double(255-IFshow2).*CATmask.*double(IVshow);
% IallCorrect(imdilate(IallCorrect>30000,ones(5))) = 0;
% IallCorrect = IallCorrect.*(1-imdilate(lungMask,ones(5))).*(1-imdilate(voidsMask,ones(5)));
% IallCorrect = IallCorrect.*(1-EATmask);
% if showProgress && ~batchShow
%     figure, imshow(IallCorrect,[]);
% end
% EATmax = max(IallCorrect(:));
% EATthreshold = 0.25*EATmax;
% IallCorrect = double(255-IFshow2).*CATmask.*double(IVshow);
% IallCorrect(IallCorrect<EATthreshold) = 0;
% IallCorrect = logical(IallCorrect);
% %IallCorrect(IVshow>0.33*EATmax) = CATmask(IVshow>0.33*EATmax);
% EATmask = logical(imdilate(EATmask,ones(3)) + IallCorrect + heartMask);
EATmask = logical(imdilate(EATmask,SEdiamond5x5) + imdilate(heartMask,SEdiamond5x5));
%figure, imshow(EATmask,[])
EATmask = imclose(EATmask,ones(7));
%figure, imshow(EATmask,[])
EATmask = FindLargestArea(EATmask);
EATmask = EATmask.*CATmask;

% fill in any holes in EAT and heart mask
EATheart = imfill(imclose(logical(EATmask+heartMask),SEdiamond5x5),'holes');
heartMask = logical(heartMask + organsMask.*EATheart);
EATmask = logical(EATmask + VATmask.*EATheart + CATmask.*EATheart);
EATmask = EATmask.*logical(1-PAATmask);

EAT = IFshow2.*uint8(EATmask);
if showProgress && ~batchShow
    figure, imshow(EAT);
end

% correct CAT mask
CATmask = CATmask.*logical(1-EATmask).*logical(1-PAATmask);
CAT   = IFshow2.*uint8(CATmask);
if showProgress && ~batchShow
    figure, imshow(CAT);
end

% correct VAT and organ mask
VATmask = VATmask.*logical(1-CATmask).*logical(1-PAATmask).*logical(1-EATmask);
organsMask = organsMask.*logical(1-heartMask).*logical(1-aortaMask);

VAT = uint8(VATmask).*IFshow2;
organs = IWshow4.*uint8(organsMask);


%% segment IMAT and VAT 
IMATmask = VATmask;

% remove spinal AT
xS = [max(1,round(size(innerMask,2)/3.0)-i) round(size(innerMask,2)/2.0) min(size(innerMask,2),round(2.0*size(innerMask,2)/3.0)+i) max(1,round(size(innerMask,2)/3.0)-i)];
yS = [size(innerMask,1) round(size(innerMask,1)/2.0) size(innerMask,1) size(innerMask,1)];
IallCorrect = imdilate(poly2mask(xS,yS,size(innerMask,1),size(innerMask,2)),SEdiamond13x13);

if showProgress && ~batchShow
    figure, imshow(IallCorrect,[]);
end

VATmask = VATmask.*(1-IallCorrect).*imerode(innerMask,ones(13));
VATmask(PAATposterior+21:end,:) = false; % removes spinal AT
VATmask = imopen(FindLargestArea(VATmask,'multiple',voidThreshold/40.0),ones(3));  %/80.0
if showProgress && ~batchShow
    figure, imshow(VATmask);
end

% reduce to VAT enclosed area
%IallCorrect = imopen(FindLargestArea(imerode((VATmask+voidsMask+lungMask).*imerode(innerMask,ones(7)),SEdiamond3x3),'multiple',voidThreshold/80.0),ones(3));   
%VATmask = imclose(FindLargestArea(imdilate(IallCorrect,ones(7)).*(VATmask+voidsMask+lungMask),'multiple',voidThreshold/80.0),ones(3));     
%VATmask(PAATposterior+21:end,:) = false; % removes spinal AT

%IallCorrect = VATmask.*IMATmask;
%if sum(lungMask(:)) > sum(IallCorrect(:))
    % VAT enclosed area defered to lungs
%    IallCorrect = imopen(FindLargestArea(imerode((CATmask+EATmask+voidsMask+lungMask+heartMask+aortaMask+PAATmask).*imerode(innerMask,ones(7)),ones(3)),'multiple',voidThreshold/80.0),ones(3));   
%    VATmask = imclose(FindLargestArea(imdilate(IallCorrect,ones(9)).*(CATmask+EATmask+voidsMask+imdilate(lungMask,ones(3))+ heartMask + aortaMask + PAATmask),'multiple',voidThreshold/80.0),ones(3));     
    %VATmask(PAATposterior+1:end,:) = false;
    %if showProgress && ~batchShow
%        figure, imshow(VATmask,[]);
    %end
%end
IallCorrect = VATmask.*IMATmask;
if (sum(lungMask(:)) > sum(IallCorrect(:))) && (sum(lungMask(:)) >  voidThreshold)
    % VAT enclosed area defered to lungs    
    IallCorrect = imopen(FindLargestArea(imerode((voidsMask+lungMask+aortaMask+PAATmask+heartMask+CATmask+EATmask).*imerode(innerMask,ones(7)),SEdiamond3x3),'multiple',voidThreshold/80.0),ones(3));   
else
    IallCorrect = imopen(FindLargestArea(imerode((VATmask+voidsMask+lungMask+aortaMask+PAATmask+heartMask+CATmask+EATmask).*imerode(innerMask,ones(7)),SEdiamond3x3),'multiple',voidThreshold/80.0),ones(3));   
end
VATfill = imfill(imdilate(IallCorrect,ones(9)),'holes');
VATmask = imclose(FindLargestArea(VATfill.*(VATmask+voidsMask+lungMask+aortaMask+PAATmask+heartMask+CATmask+EATmask),'multiple',voidThreshold/80.0),ones(3));     
if showProgress && ~batchShow
    figure, imshow(VATmask);
end

% upper convex hull mask 
IallCorrect(:) = false;
IallCorrect(1:aortaAnterior,:) = true;

if showProgress && ~batchShow
    figure, imshow(IallCorrect,[]);
end

statsVAT = regionprops(VATmask.*IallCorrect, 'ConvexHull');
if size(statsVAT,1)>0
    VATfill = logical(VATfill + imopen(roipoly(VATmask.*IallCorrect,statsVAT(1).ConvexHull(:,1),statsVAT(1).ConvexHull(:,2)),ones(3)));
    VATmask = logical(VATmask + VATfill.*IMATmask);
    VATmask = FindLargestArea(VATmask,'multiple',voidThreshold/40.0); %/80.0
    if showProgress && ~batchShow
        figure, imshow(VATmask);
    end
end

% left lobe convex hull mask
xL = [round(size(innerMask,2)/3.0) round(size(innerMask,2)/2.0) round(size(innerMask,2)/3.0) 1 1 round(size(innerMask,2)/3.0)];
yL = [1 aortaAnterior size(innerMask,1) size(innerMask,1) 1 1];
IallCorrect = poly2mask(xL,yL,size(innerMask,1),size(innerMask,2));
if showProgress && ~batchShow
    figure, imshow(IallCorrect,[]);
end

statsVAT = regionprops(VATmask.*IallCorrect, 'ConvexHull');
if size(statsVAT,1)>0
VATfill = logical(VATfill + imdilate(roipoly(VATmask.*IallCorrect,statsVAT(1).ConvexHull(:,1),statsVAT(1).ConvexHull(:,2)),ones(3)));
VATmask = logical(VATmask + VATfill.*IMATmask);
    if showProgress && ~batchShow
       figure, imshow(VATmask);
    end
end

% right lobe convex hull mask
xR = [round(2.0*size(innerMask,2)/3.0) round(size(innerMask,2)/2.0) round(2.0*size(innerMask,2)/3.0) size(innerMask,2) size(innerMask,2) round(2.0*size(innerMask,2)/3.0)];
yR = [1 aortaAnterior size(innerMask,1) size(innerMask,1) 1 1];
IallCorrect = poly2mask(xR,yR,size(innerMask,1),size(innerMask,2));

if showProgress && ~batchShow
    figure, imshow(IallCorrect,[]);
end

statsVAT = regionprops(VATmask.*IallCorrect, 'ConvexHull');
if size(statsVAT,1) > 0
VATfill = logical(VATfill + imdilate(roipoly(VATmask.*IallCorrect,statsVAT(1).ConvexHull(:,1),statsVAT(1).ConvexHull(:,2)),SEdiamond3x3));   
VATmask = logical(VATmask + VATfill.*IMATmask);
if showProgress && ~batchShow
    figure, imshow(VATmask);
end
end

%VATfill = imclose(logical(VATmask+voidsMask+lungMask+heartMask+CATmask+aortaMask+PAATmask),SEdiamond3x3);
%VATfill = imclose(logical(VATmask+voidsMask+lungMask+heartMask+CATmask+aortaMask+PAATmask),SEdiamond9x9);
VATfill = imfill(imclose(VATfill,SEdiamond9x9),'holes');
VATfill = FindLargestArea(VATfill,'multiple',voidThreshold/80.0);
% remove spinal AT
xS = [max(1,round(size(innerMask,2)/3.0)-i) round(size(innerMask,2)/2.0) min(size(innerMask,2),round(2.0*size(innerMask,2)/3.0)+i) max(1,round(size(innerMask,2)/3.0)-i)];
yS = [size(innerMask,1) round(size(innerMask,1)/2.0) size(innerMask,1) size(innerMask,1)];
IallCorrect = imdilate(poly2mask(xS,yS,size(innerMask,1),size(innerMask,2)),SEdiamond13x13);
if showProgress && ~batchShow
    figure, imshow(IallCorrect,[]);
end
VATfill = VATfill.*(1-IallCorrect).*imerode(innerMask,ones(11));
%VATfill(PAATposterior+21:end,:) = false; % removes spinal AT
if i > first_i
    VATfill = logical(VATfill + imdilate((VAT3d(:,:,i-1)>0),SEdiamond3x3).*IMATmask.*imerode(innerMask,ones(13))); 
end
if showProgress && ~batchShow
    figure, imshow(VATfill);
end

% use the VATfill as the seed to grow the VAT mask
% VATmask = VATfill;
% for n = 1:10
%     VATmask = imdilate(VATmask,SEdiamond3x3);
     VATmask = VATfill.*IMATmask;
% end

VATfill = imerode(imfill(imdilate(VATmask,ones(11)),'holes'),ones(9)).*imerode(innerMask,SEdiamond9x9);

IallCorrect = IMATmask;
IMATmask = logical(IMATmask.*imerode((1-VATmask),SEdiamond3x3));
VATmask = logical(imdilate(VATmask,SEdiamond3x3).*(1-IMATmask)).*IallCorrect;
%IMATmask = logical(IMATmask.*(1-VATmask));

if showProgress && ~batchShow
    figure, imshow(IMATmask);
end

% define thoracic fat mask if necessary
if i < heartMax
   % assume all fat in body cavity is TAT since above the diaphram 
   TATmask = VATmask;
   VATmask = 0.*TATmask;
else
    TATmask = 0.*VATmask;
    if i > first_i
        % allow for TAT to continue from more 
        TATmask = imdilate((TAT3d(:,:,i-1)>0),SEdiamond3x3).*logical(VATmask + IMATmask).*imerode(innerMask,ones(13)); 
        IallCorrect = imdilate(IMATmask,ones(13)).*VATmask;
        TATmask = logical(TATmask + IallCorrect);
        % mask with lung region
        IallCorrect = (lungMask+heartMask+CATmask+EATmask);
        if sum(sum(IallCorrect)) > 200 
            % do nothing
        else
            TATmask = 0.*VATmask;
        end
    end
    VATmask = VATmask.*logical(1-TATmask);
    IMATmask = IMATmask.*logical(1-TATmask);
end
TAT  = uint8(TATmask).*IFshow2;
VAT  = uint8(VATmask).*IFshow2;
IMAT = uint8(IMATmask).*IFshow2;

%% tally area sums

SCATarea  = sum(sum(SCATmask))
TATarea   = sum(sum(TATmask))
VATarea   = sum(sum(VATmask))
IMATarea  = sum(sum(IMATmask))
organArea = sum(sum(organsMask))
voidArea  = sum(sum(voidsMask))

SCATfuzzyArea  = sum(sum(SCAT))/255.
TATfuzzyArea   = sum(sum(TAT))/255.
VATfuzzyArea   = sum(sum(VAT))/255.
IMATfuzzyArea  = sum(sum(IMAT))/255.
organFuzzyArea = sum(sum(organs))/255.
voidFuzzyArea  = sum(sum(voids))/255.

% thoracic tissue array
heart = IWshow4.*uint8(heartMask);

tissues(:,:,1) = 0.667*organs + SCAT + aorta +  heart + EAT + 0.667*CAT;
tissues(:,:,2) = TAT + VAT + SCAT + IMAT + PAAT + 0.125*heart + CAT + EAT;
tissues(:,:,3) = 0.5*voids + aorta + PAAT + 0.125*heart + lungs + EAT + 0.333*CAT;

if ~batchShow
figure, imshow(tissues);
end

tissues(:,:,1) = 0.667*organs + SCAT + IMAT + aorta +  heart + EAT + IMAT + 0.667*CAT;
tissues(:,:,2) = 0.667*TAT + VAT + SCAT + 0.667*IMAT + PAAT + 0.125*heart + CAT + EAT;
tissues(:,:,3) = 0.5*voids + 0.333*IMAT + aorta + PAAT + 0.125*heart + lungs + EAT + 0.333*CAT;

if i == slice2show
figure, imshow(tissues);
end

% lung-heart-fat array
TTarray = [lungs,heart,CAT];
%if showProgress && ~batchShow
%figure, montage(TTarray);
%end

TTtissues(:,:,1) = heart + aorta + EAT + 0.667*CAT;
TTtissues(:,:,2) = 0.125*heart + TAT + VAT + IMAT + CAT + PAAT + EAT;
TTtissues(:,:,3) = 0.125*heart + lungs + aorta + PAAT + EAT + 0.333*CAT;

if ~batchShow
figure, imshow(TTtissues);
end

TTtissues(:,:,1) = heart + aorta + EAT + IMAT + 0.667*CAT;
TTtissues(:,:,2) = 0.125*heart + VAT + CAT + PAAT + EAT +0.667*IMAT + 0.667*TAT;
TTtissues(:,:,3) = 0.125*heart + lungs + aorta + PAAT + EAT + 0.333*CAT + 0.333*IMAT;

if ~batchShow
figure, imshow(TTtissues);
end

TTtissues(:,:,1) = heart + aorta + EAT + 0.667*CAT;
TTtissues(:,:,2) = 0.125*heart + CAT + PAAT + EAT;
TTtissues(:,:,3) = 0.125*heart + lungs + aorta + PAAT + EAT + 0.333*CAT;

if ~batchShow
figure, imshow(TTtissues);
end

heartTissues(:,:,1) = heart + EAT + aorta + 0.667*CAT;
heartTissues(:,:,2) = 0.125*heart + CAT + PAAT + EAT;
heartTissues(:,:,3) = 0.125*heart + aorta + PAAT + EAT + 0.333*CAT;

if ~batchShow
figure, imshow(heartTissues);
end

SCATvolume(i)  = double(sum(sum(SCATmask)))*voxelVolume/1000.0;  % [mL]
TATvolume(i)   = double(sum(sum(TATmask)))*voxelVolume/1000.0;   % [mL]
VATvolume(i)   = double(sum(sum(VATmask)))*voxelVolume/1000.0;   % [mL]
IMATvolume(i)  = double(sum(sum(IMATmask)))*voxelVolume/1000.0;  % [mL]
organsVolume(i)= double(sum(sum(organsMask)))*voxelVolume/1000.0;% [mL]
voidsVolume(i) = double(sum(sum(voidsMask)))*voxelVolume/1000.0; % [mL]
lungVolume(i)  = double(sum(sum(lungMask)))*voxelVolume/1000.0;  % [mL]
heartVolume(i) = double(sum(sum(heartMask)))*voxelVolume/1000.0; % [mL]
aortaVolume(i) = double(sum(sum(aortaMask)))*voxelVolume/1000.0; % [mL]
CATvolume(i)   = double(sum(sum(CATmask)))*voxelVolume/1000.0;   % [mL]
PAATvolume(i)  = double(sum(sum(PAATmask)))*voxelVolume/1000.0;  % [mL]
EATvolume(i)   = double(sum(sum(EATmask)))*voxelVolume/1000.0;   % [mL]

SCAT3d(:,:,i)   = SCAT;
TAT3d(:,:,i)    = TAT;
VAT3d(:,:,i)    = VAT;
IMAT3d(:,:,i)   = IMAT;
organs3d(:,:,i) = organs;
voids3d(:,:,i)  = voids;
lung3d(:,:,i)   = lungs;
heart3d(:,:,i)  = heart;
aorta3d(:,:,i)  = aorta;
CAT3d(:,:,i)    = CAT;
PAAT3d(:,:,i)   = PAAT;
EAT3d(:,:,i)    = EAT;

else
PAATmask = 0*aortaMask;
PAAT = IFshow2.*uint8(PAATmask); 
%VATmask = VATmask.*logical(1-PAATmask);    
    
% % auto detect aorta
% % finds the circles
% [r c rad] = circlefinder(IWshow2,3,12,[],[],true);
% 
% % draws the circles
% for n=1:length(rad)
%     imcircle = RGBCircle(IWshow2,r(n),c(n),rad(n), [255 255 255], 2);
% end
% figure, imshow(imcircle);
imcircle = 0;
    

%% segment IMAT and VAT 

IMATmask = VATmask;

if i < FH_i
% reduce to VAT enclosed area
IallCorrect = imopen(FindLargestArea(imerode((VATmask+voidsMask).*imerode(innerMask,ones(7)),ones(3)),'multiple',voidThreshold/80.0),ones(3));   
VATmask = imclose(FindLargestArea(imdilate(IallCorrect,ones(7)).*(VATmask+voidsMask),'multiple',voidThreshold/80.0),ones(3));     

if showProgress && ~batchShow
    figure, imshow(VATmask,[]);
end

% upper convex hull mask 
IallCorrect(:) = false;
IallCorrect(1:aortaAnterior,:) = true;

if showProgress && ~batchShow
    figure, imshow(IallCorrect,[]);
end

statsVAT = regionprops(VATmask.*IallCorrect, 'ConvexHull');
if size(statsVAT,1) > 0
VATmask = logical(VATmask + imopen(roipoly(VATmask.*IallCorrect,statsVAT(1).ConvexHull(:,1),statsVAT(1).ConvexHull(:,2)).*IMATmask,ones(3)));
VATmask = FindLargestArea(VATmask,'multiple',voidThreshold/80.0);
end

VATanterior = find(sum(VATmask,2)>0,1);
IallCorrect(1:VATanterior,:) = false;
statsVAT = regionprops(VATmask.*IallCorrect, 'ConvexHull');
if size(statsVAT,1) > 0
VATmask = logical(VATmask + imdilate(roipoly(VATmask.*IallCorrect,statsVAT(1).ConvexHull(:,1),statsVAT(1).ConvexHull(:,2)).*IMATmask,SEdiamond3x3));
end
if showProgress && ~batchShow
    figure, imshow(VATmask);
end

% left lobe convex hull mask
xL = [round(size(innerMask,2)/3.0) round(size(innerMask,2)/2.0) round(size(innerMask,2)/3.0) 1 1 round(size(innerMask,2)/3.0)];
yL = [1 aortaAnterior size(innerMask,1) size(innerMask,1) 1 1];
IallCorrect = poly2mask(xL,yL,size(innerMask,1),size(innerMask,2));

if showProgress && ~batchShow
    figure, imshow(IallCorrect,[]);
end

statsVAT = regionprops(VATmask.*IallCorrect, 'ConvexHull');
if size(statsVAT,1) > 0
VATmask = logical(VATmask + imdilate(roipoly(VATmask.*IallCorrect,statsVAT(1).ConvexHull(:,1),statsVAT(1).ConvexHull(:,2)),ones(3)).*IMATmask);
end
if showProgress && ~batchShow
    figure, imshow(VATmask);
end

% right lobe convex hull mask
xR = [round(2.0*size(innerMask,2)/3.0) round(size(innerMask,2)/2.0) round(2.0*size(innerMask,2)/3.0) size(innerMask,2) size(innerMask,2) round(2.0*size(innerMask,2)/3.0)];
yR = [1 aortaAnterior size(innerMask,1) size(innerMask,1) 1 1];
IallCorrect = poly2mask(xR,yR,size(innerMask,1),size(innerMask,2));

if showProgress && ~batchShow
    figure, imshow(IallCorrect,[]);
end

statsVAT = regionprops(VATmask.*IallCorrect, 'ConvexHull');
if size(statsVAT,1) > 0
    VATmask = logical(VATmask + imdilate(roipoly(VATmask.*IallCorrect,statsVAT(1).ConvexHull(:,1),statsVAT(1).ConvexHull(:,2)),ones(3)).*IMATmask).*imerode(innerMask,SEdiamond3x3).*logical(1-voidsMask);
end

if showProgress && ~batchShow
    figure, imshow(VATmask);
end

% remove spinal AT
if i < FH_i
xS = [max(1,round(size(innerMask,2)/3.0)-i) round(size(innerMask,2)/2.0) min(size(innerMask,2),round(2.0*size(innerMask,2)/3.0)+i) max(1,round(size(innerMask,2)/3.0)-i)];
yS = [size(innerMask,1) round(size(innerMask,1)/2.0) size(innerMask,1) size(innerMask,1)];
IallCorrect = imdilate(poly2mask(xS,yS,size(innerMask,1),size(innerMask,2)),SEdiamond13x13);
if i > (FH_i - 50)
    IallCorrect = imerode(IallCorrect,ones(2*(i-(FH_i - 50))));
end
if showProgress && ~batchShow
    figure, imshow(IallCorrect,[]);
end
VATmask = VATmask.*(1-IallCorrect);
else
    % do nothing
end

if showProgress && ~batchShow
    figure, imshow(VATmask);
end

if i > FH_i-30
    PCwidth = 0.2*(FH_i+8-i)^2;
    PCmask = VATmask;
    PCmask(:) = true;
    PCmask(:,1:round(0.5*size(PCmask,2)-PCwidth)) = false;
    PCmask(:,round(0.5*size(PCmask,2)+PCwidth):end) = false;
    VATmask = VATmask.*PCmask;
    clear 'PCmask';
end

else
    VATmask = 0.*IMATmask;
end

if (i > first_i) && (i < FH_i + 3)
    if i < FH_i-30
        VATmask = logical(VATmask + (VAT3d(:,:,i-1)>0).*IMATmask.*imerode(innerMask,ones(11)));
    else
        VATmask = logical(VATmask + imopen((VAT3d(:,:,i-1)>0),SEdiamond5x5).*imdilate(VATmask,SEdiamond7x7).*imerode(innerMask,ones(11)));
    end
end
VATmask = VATmask.*imerode(innerMask,SEdiamond7x7);

VATfill = imerode(imfill(imdilate(VATmask,ones(11)),'holes'),ones(9));

IallCorrect = IMATmask;
IMATmask = logical(IMATmask.*imerode((1-VATmask),SEdiamond3x3));
VATmask = logical(imdilate(VATmask,SEdiamond3x3).*(1-IMATmask)).*IallCorrect;

if showProgress && ~batchShow
    figure, imshow(IMATmask);
end

%% associate unassigned pixels

% needed for the remaining pixels not assigned to any mask
voidsCorrect = logical(bwmorph(bodyMask,'erode',3) - logical(SCATmask + VATmask + IMATmask + organsMask +voidsMask));
% intensities of voxels near voids (fat saturated signal) 
NvoidsF = IFshow2.*uint8(voidsCorrect);
% intensities of voxels near voids (water saturated signal) 
NvoidsW = IWshow4.*uint8(voidsCorrect);
% intensities of voxels near voids (void saturated signal) 
NvoidsV = voids.*uint8(voidsCorrect);
% unassigned pixels are assigned based on the higher signal
organsCorrect = ((NvoidsW) > (NvoidsF))&((NvoidsW) > (NvoidsV)); 
fatCorrect = ((NvoidsF) > (NvoidsW))&((NvoidsF) > (NvoidsV)); 
voidsCorrect = logical(voidsCorrect - organsCorrect - fatCorrect);
% add unassigned voxels to appropriate masks
VATcorrect = logical((fatCorrect).*VATfill);
IMATcorrect = logical((fatCorrect).*innerMask.*(1 - VATfill));
SCATcorrect = logical(fatCorrect.*(1 - VATcorrect).*(1 - IMATcorrect));
% update SCAT, VAT and organ masks
SCATmask(SCATcorrect) = 1;
SCAT = uint8(SCATmask).*IFshow2;
VATmask(VATcorrect) = 1;
VAT = uint8(VATmask).*IFshow2;
IMAT = uint8(IMATmask).*IFshow2;
organsMask(organsCorrect) = 1;
organs = IWshow4.*uint8(organsMask);
voidsMask(organsCorrect) = 0;
voidsMask(fatCorrect) = 0;
voidsMask(voidsCorrect) = 1;
voids = uint8((peakPixelBinW/255.0)*(255-IWshow4).*uint8(voidsMask)); 

organsMask = organsMask.*(1-aortaMask);

% adipose tissue array
ATarray = [SCAT,VAT,IMAT,organs,voids];

tissues(:,:,1) = 0.667*organs + SCAT + aorta;
tissues(:,:,2) = VAT +IMAT + SCAT + PAAT;
tissues(:,:,3) = 0.5*voids + aorta +PAAT;

if ~batchShow
figure, imshow(tissues);
end

tissues(:,:,1) = 0.667*organs + SCAT + IMAT + aorta;
tissues(:,:,2) = VAT + SCAT + 0.667*IMAT + PAAT;
tissues(:,:,3) = 0.5*voids + 0.333*IMAT + aorta +PAAT;

if i == slice2show
figure, imshow(tissues);
end

% if i == (length(upperSlices) + 1)
% 
% SCATvolume(i-1)  = 0.5*(SCATvolume(i-1) + double(sum(sum(SCATmask)))*voxelVolume/1000.0);  % [mL]
% VATvolume(i-1)   = 0.5*(VATvolume(i-1)  + double(sum(sum(VATmask)))*voxelVolume/1000.0);   % [mL]
% organsVolume(i-1)= 0.5*(organsVolume(i-1)+double(sum(sum(organsMask)))*voxelVolume/1000.0);% [mL]
% voidsVolume(i-1) = 0.5*(voidsVolume(i-1) +double(sum(sum(voidsMask)))*voxelVolume/1000.0); % [mL]
% 
% SCAT3d(:,:,i-1)  = 0.5*(SCAT3d(:,:,i-1)+SCAT);
% VAT3d(:,:,i-1)   = 0.5*(VAT3d(:,:,i-1)+VAT);
% organs3d(:,:,i-1)= 0.5*(organs3d(:,:,i-1)+organs);
% voids3d(:,:,i-1) = 0.5*(voids3d(:,:,i-1)+voids);
% 
% elseif i > (length(upperSlices) + 1)
%     
% SCATvolume(i-1)  = double(sum(sum(SCATmask)))*voxelVolume/1000.0;  % [mL]
% VATvolume(i-1)   = double(sum(sum(VATmask)))*voxelVolume/1000.0;   % [mL]
% organsVolume(i-1)= double(sum(sum(organsMask)))*voxelVolume/1000.0;% [mL]
% voidsVolume(i-1) = double(sum(sum(voidsMask)))*voxelVolume/1000.0; % [mL]
% 
% SCAT3d(:,:,i-1)  = SCAT;
% VAT3d(:,:,i-1)   = VAT;
% organs3d(:,:,i-1)= organs;
% voids3d(:,:,i-1) = voids;
% 
% else

SCATvolume(i)  = double(sum(sum(SCATmask)))*voxelVolume/1000.0;  % [mL]
VATvolume(i)   = double(sum(sum(VATmask)))*voxelVolume/1000.0;   % [mL]
IMATvolume(i)  = double(sum(sum(IMATmask)))*voxelVolume/1000.0;   % [mL]
organsVolume(i)= double(sum(sum(organsMask)))*voxelVolume/1000.0;% [mL]
voidsVolume(i) = double(sum(sum(voidsMask)))*voxelVolume/1000.0; % [mL]

SCAT3d(:,:,i)  = SCAT;
VAT3d(:,:,i)   = VAT;
IMAT3d(:,:,i)  = IMAT;
organs3d(:,:,i)= organs;
voids3d(:,:,i) = voids;    
    
%end
end

end

% properly orient resulting saggital image
Isaggital = rot90(Isaggital',2); 

varlist = {'niiFup_reslice', 'niiFlo_reslice', 'niiWup_reslice', 'niiWlo_reslice','aortaSeed','IvetoCorrect','statsAorta','statsHeart','statsEAT','statsInner',...
           'SCATmask','VATmask','TATmask','organsMask','voidsMask','bodyMask','innerMask','heartMask','CATmask','PAATmask','aortaMask','CATheart','SCATinnerMask','leftLungVetoMask',...
           'cutMask','lungMask','SCATplusExternalMask','nonATmask','SCATfill','perimeterMask','externalMask','centroidMask','linkMask','skinMask','EATmask','IMATmask',...
           'SCAT','VAT','organs','voids','heart','lungs','CAT','EAT','PAAT','ATarray','TTarray','TTtissues','heartTissues','aorta','notLung','nonAT','tissues','IMAT','IMATcorrect',...
           'voidsCorrect','organsCorrect','fatCorrect','VATcorrect','SCATcorrect','PAATcorrect','heartCorrect','diaphramCorrect','SCATclosed','SCATdilate','SCATinnerMaskFill',...
           'NvoidsF','NvoidsW','NvoidsV','CATveto','lungVeto','aortaVeto','cutSaveLine','seg_IFbin','seg_IWbin','erodedBWfat','filledBWfat','openedBWfat','labeledBWaorta','imcircle'...
           'heartBase','aortaBase','PAATbase','heartFatBase','CATbase','voidsBase','fillbase','SCATbase','SCATbase2','filledFatBase','innerBase','lungBase','EATheart',...
           'IFraw','IFshow','IFshow1','IFshow2','IFcorrect','IWraw','IWshow','IWshow2','IWshow3','IWshow4','IWcorrect','Iplus','IV','IVshow','IFrawVeto','IWrawVeto','IWrawSample','IFrawSample','IVraw'};
% 
if ~showProgress && writeResults       
    clear(varlist{:})
    clear('varlist')
    save(resultsFile);
end