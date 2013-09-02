function threshed_image = find_threshed_image(high_passed_image, filter_thresh, ...
    proximity_filter,min_seed_size)

if (length(filter_thresh) == 1)
    threshed_image = high_passed_image >= filter_thresh;
else
    high_threshed_image = high_passed_image >= filter_thresh(2);
    
    if (not(isnan(min_seed_size)))
        high_threshed_labels = bwlabel(high_threshed_image);
        high_threshed_props = regionprops(high_threshed_labels,'Area'); %#ok<MRPBW>
        
        high_threshed_image = ismember(high_threshed_labels,find([high_threshed_props.Area] >= 4));
    end
    
    high_threshed_image = imdilate(high_threshed_image,strel('disk',proximity_filter));
    
    low_threshed_image = high_passed_image >= filter_thresh(1);
    low_thresh_bwlabel = bwlabel(low_threshed_image,4);
    
    overlap_labels = unique(low_thresh_bwlabel.*high_threshed_image);
    if (overlap_labels(1) == 0)
        overlap_labels = overlap_labels(2:end);
    end
    
    threshed_image = ismember(low_thresh_bwlabel,overlap_labels);
end
end
