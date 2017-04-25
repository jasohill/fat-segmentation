function valleyPixelBin = AutorestrictThreshold(I,valleyPixelBin,nbrThresholds)
% restrict foreground threshold to be equal to or less than the next to
% highest automatically determined threshold via Otsu's method
thresh = multithresh(I,nbrThresholds);

% restrict foreground threshold to be equal to or less than the next to
% highest automatically determined threshold
if valleyPixelBin > thresh(nbrThresholds-1)
    valleyPixelBin = thresh(nbrThresholds-1);
end
if valleyPixelBin < thresh(2)
    valleyPixelBin = thresh(2);
end

end

