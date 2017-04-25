function [peakPixelBin,valleyPixelBin] = FindHistogramExtrema(I,backgroundThreshold,foregroundThreshold,showFlag);
% Find the peak of the smoothed histogram of I (pixel with the smallest bin size)

% Find histogram of the image foreground
[counts,binLocations] = imhist(I(I>backgroundThreshold));
% smooth the histogram bin counts with a median filter
medCounts = medfilt1(counts,7);
% remove influence of background from the smoothed histogram 
medCounts(1:foregroundThreshold) = 0;
medCounts(foregroundThreshold+1:253) = medCounts(foregroundThreshold+2:254);
medCounts(254) = medCounts(255);
medCounts(255) = 0;

if showFlag
    figure, stem(binLocations,medCounts);
end

% find the peak of the smoothed histogram bins
maxIndeces0 = find(medCounts == max(medCounts));
maxIndex0 = ceil(median(maxIndeces0));
if maxIndex0 == 1
    maxIndex0 = 2;
end
if maxIndex0 == length(counts)
    maxIndex0 = counts-1;
end
count3 = [counts(maxIndex0-1),counts(maxIndex0),counts(maxIndex0+1)];
maxIndeces0 = find(count3 == max(count3));

peakPixelBin = maxIndex0 - 3 + ceil(median(maxIndeces0));

% fill in the smoothed histogram to enable finding the valley
medCounts(1:foregroundThreshold) = medCounts(maxIndex0);
medCounts(maxIndex0:end) = medCounts(maxIndex0);

if showFlag
    figure, stem(binLocations,medCounts);
end

% find the valley of the smoothed histogram bins
minIndeces0 = find(medCounts == min(medCounts));
minIndex0 = ceil(median(minIndeces0));  
count3 = [counts(minIndex0-1),counts(minIndex0),counts(minIndex0+1)];
minIndeces0 = find(count3 == min(count3));

valleyPixelBin = minIndex0 - 3 + median(minIndeces0);

valleyPixelBin =  round(valleyPixelBin*(255/peakPixelBin));

end