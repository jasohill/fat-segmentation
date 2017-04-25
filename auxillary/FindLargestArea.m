function BWout = FindLargestArea(BWin,varargin)
% label separated connected regions and return one with largest area

% handle variable inputs
if nargin == 3
    code = varargin{1};
    threshold = varargin{2};
else
    code = [];
end

[labeledBW,nbr_labels] = bwlabel(BWin,4);

stats = regionprops(labeledBW,'area');
% identify the connected region with the largest area
label_areas = cat(1, stats.Area);

if isempty(code)
    k_largest = max(label_areas);
    
    label = find(label_areas == k_largest);
    
    if isempty(label)
        BWout = 0.*BWin;
    else
        if length(label) > 1
            BWout = 0.*BWin;            
            for l=1:length(label)
                BWout = logical(BWout + (labeledBW == label(l))); 
            end
                
        else
            BWout = (labeledBW == label);
        end
    end
    
elseif strcmp(code,'multiple')
    labels = find(label_areas > threshold);

    BWout = 0.*BWin;
    if isempty(labels)
        % do nothing
    else
        % build output
        for label = 1:length(labels) 
            BWout = BWout|(labeledBW == labels(label));
        end
    end
end



end

