clear all; close all;

filenameBase = 'subject02_f';
filenameIn   = [filenameBase,'.mat'];
filenameOut  = [filenameBase,'.xlsx'];

load(filenameIn)

i = 1:Nslices;
i(Nslices+1:Nslices+2)      = NaN;
levels(Nslices+1:Nslices+2) = NaN;
slices(Nslices+1:Nslices+2) = NaN;
voxelVolumes(Nslices+1:Nslices+2) = NaN;
SCATvolume(Nslices+1)   = sum(SCATvolume(diaphram:Nslices));
SCATvolume(Nslices+2)   = sum(SCATvolume(1:Nslices));
VATvolume(Nslices+1)    = sum(VATvolume(diaphram:Nslices));
VATvolume(Nslices+2)    = sum(VATvolume(1:Nslices));
organsVolume(Nslices+1) = sum(organsVolume(diaphram:Nslices));
organsVolume(Nslices+2) = sum(organsVolume(1:Nslices));
voidsVolume(Nslices+1)  = sum(voidsVolume(diaphram:Nslices));
voidsVolume(Nslices+2)  = sum(voidsVolume(1:Nslices));
lungVolume(Nslices+1)   = sum(lungVolume(diaphram:Nslices));
lungVolume(Nslices+2)   = sum(lungVolume(1:Nslices));
heartVolume(Nslices+1)  = sum(heartVolume(diaphram:Nslices));
heartVolume(Nslices+2)  = sum(heartVolume(1:Nslices));
aortaVolume(Nslices+1)  = sum(aortaVolume(diaphram:Nslices));
aortaVolume(Nslices+2)  = sum(aortaVolume(1:Nslices));
CATvolume(Nslices+1)    = sum(CATvolume(diaphram:Nslices));
CATvolume(Nslices+2)    = sum(CATvolume(1:Nslices));
PAATvolume(Nslices+1)   = sum(PAATvolume(diaphram:Nslices));
PAATvolume(Nslices+2)   = sum(PAATvolume(1:Nslices));

% Create volume count table by slice [cc]
VolumeCounts = table(i',levels',slices',voxelVolumes',SCATvolume',VATvolume',...
                    organsVolume',voidsVolume',lungVolume',heartVolume',...
                    aortaVolume',CATvolume',PAATvolume',...
                    'VariableNames',{'CombinedSlice','Scan',...
                    'OriginalSlice','VoxelVolume','SCAT','VAT',...
                    'organs','voids','lungs','heart','aorta','CAT','PAAT'});

% Create normalized intensity sum count table by slice (that is winner take all) 
for n = 1:Nslices-1
    SCATintensity(n)   = double(sum(sum(SCAT3d(:,:,n))))/255.0;
    VATintensity(n)    = double(sum(sum(VAT3d(:,:,n))))/255.0;
    organsIntensity(n) = double(sum(sum(organs3d(:,:,n))))/255.0;
    voidsIntensity(n)  = double(sum(sum(voids3d(:,:,n))))/255.0;    
    if lungVolume(n) > 0
        lungIntensity(n) = double(sum(sum(lung3d(:,:,n))))/255.0;  
    else
        lungIntensity(n) = 0;
    end
    if heartVolume(n) > 0
        heartIntensity(n) = double(sum(sum(heart3d(:,:,n))))/255.0;  
    else
        heartIntensity(n) = 0;
    end
    if aortaVolume(n) > 0
        aortaIntensity(n) = double(sum(sum(aorta3d(:,:,n))))/255.0;  
    else
        aortaIntensity(n) = 0;
    end
    if CATvolume(n) > 0
        CATintensity(n) = double(sum(sum(CAT3d(:,:,n))))/255.0;  
    else
        CATintensity(n) = 0;
    end
    if PAATvolume(n) > 0
        PAATintensity(n) = double(sum(sum(PAAT3d(:,:,n))))/255.0;  
    else
        PAATintensity(n) = 0;
    end       
end
SCATintensity(Nslices)   = sum(SCATintensity);
VATintensity(Nslices)    = sum(VATintensity);
organsIntensity(Nslices) = sum(organsIntensity);
voidsIntensity(Nslices)  = sum(voidsIntensity);
lungIntensity(Nslices)   = sum(lungIntensity);
heartIntensity(Nslices)  = sum(heartIntensity);
aortaIntensity(Nslices)  = sum(aortaIntensity);
CATintensity(Nslices)    = sum(CATintensity);
PAATintensity(Nslices)   = sum(PAATintensity);

n = 1:Nslices-1;
n(Nslices) = NaN;
IntensityCounts = table(n',SCATintensity',VATintensity',...
                    organsIntensity',voidsIntensity',lungIntensity',heartIntensity',...
                    aortaIntensity',CATintensity',PAATintensity',...
                    'VariableNames',{'sliceNumber','SCAT','VAT',...
                    'organs','voids','lungs','heart','aorta','CAT','PAAT'});

writetable(IntensityCounts,filenameOut,'Sheet',2);
writetable(VolumeCounts,filenameOut,'Sheet',1);