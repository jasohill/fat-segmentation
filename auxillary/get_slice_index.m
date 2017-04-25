function slice = get_slice_index(slices,levels,index,level)
% The local function get_slice_index returns the corresponding slice index to i in level.
    slice = [];
    
    if strcmp(level,'upper')
        level = 1;
    elseif strcmp(level,'lower')
        level = 2;
    end
    if strcmp(level,'top')
        level = 1;
    elseif strcmp(level,'bottom')
        level = 2;
    end
    indeces = find(slices == index);

    if isempty(indeces)
        warning('Requested slice index not found!')
    elseif length(indeces) == 1
        if levels(indeces) == level
            slice = indeces;
        else
            warning('Requested slice index not found!')
        end
    else
        if levels(indeces(1)) == level
            slice = indeces(1);
        elseif levels(indeces(2)) == level
            slice = indeces(2);
        else
            warning('Requested slice index not found!')            
        end
    end
end