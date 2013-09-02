test_image = double(imread('../data/kiehart/mov1/08k21_30000.tif'));

%High-pass filtering
I_filt = fspecial('disk',10);
blurred_image = imfilter(test_image,I_filt,'replicate');
high_passed_image = test_image - blurred_image;

%Thresholding/Largest Object Selection
threshed_image = find_threshed_image(high_passed_image, [.1,20],0,1);
regions = bwlabel(threshed_image);
props = regionprops(regions,'Area'); %#ok<MRPBW>

cell_outlines = ismember(regions,find([props.Area] == max([props.Area])));
cell_outlines = bwmorph(cell_outlines,'skel',Inf);

cell_bodies_binary = not(cell_outlines); 
cell_bodies = bwlabel(cell_bodies_binary,4);

%Filter objects/cell bodies at edge
edge_pixels = [cell_bodies(1,:),cell_bodies(end,:), ...
    cell_bodies(:,1)',cell_bodies(:,end)'];
edge_labels = unique(edge_pixels);
edge_labels = edge_labels(edge_labels ~= 0);
not_edge_labels = setdiff(unique(cell_bodies(:)),edge_labels);
not_edge_labels = not_edge_labels(not_edge_labels ~= 0);

%Filter small objects
cell_bodies = ismember(cell_bodies,not_edge_labels).*cell_bodies;
cell_props = regionprops(cell_bodies,'Area');
cell_bodies = ismember(cell_bodies,find([cell_props.Area] > 50)).*cell_bodies;


% Renumber adhesions to be sequential
cell_nums = unique(cell_bodies);
assert(cell_nums(1) == 0, 'Background pixels not found after building adhesion label matrix')
for i = 2:length(cell_nums)
    cell_bodies(cell_bodies == cell_nums(i)) = i - 1;
end

test_image_norm = (test_image - min(test_image(:)))/range(test_image(:));
imshow(create_highlighted_image(test_image_norm,cell_bodies,'mix_percent',0.5));