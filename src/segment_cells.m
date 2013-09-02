function segment_cells(image_folder,varargin)
% SEGMENT_CELLS    Segments cell bodies from N-cadherin images.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
i_p = inputParser;
i_p.StructExpand = true;
i_p.addRequired('image_folder',@(x)exist(x,'dir') == 7);

i_p.parse(image_folder,varargin{:});

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_filenames = dir(image_folder);
image_filenames = image_filenames(3:end);

for i_num = 1:length(image_filenames)
    image = double(imread(fullfile(image_folder,image_filenames(i_num).name)));
    
    %High-pass filtering
    I_filt = fspecial('disk',10);
    blurred_image = imfilter(image,I_filt,'replicate');
    high_passed_image = image - blurred_image;
    
    %Thresholding/Largest Object Selection
    threshed_image = find_threshed_image(high_passed_image,[0,20],0,1);
    regions = bwlabel(threshed_image);
    props = regionprops(regions,'Area'); %#ok<MRPBW>
    
    cell_outlines = ismember(regions,find([props.Area] == max([props.Area])));
    cell_outlines = bwmorph(cell_outlines,'skel',Inf);
    cell_outlines = imdilate(cell_outlines,strel('disk',4));
    cell_outlines = bwmorph(cell_outlines,'skel',Inf);
    
    after = bwmorph(cell_outlines,'spur');
    while (not(all(all(cell_outlines == after))))
        cell_outlines = after;
        after = bwmorph(cell_outlines,'spur');
    end
    
    cell_bodies_binary = not(cell_outlines);
    cell_bodies = bwlabel(cell_bodies_binary,4);
    
%     %Filter objects/cell bodies at edge
%     edge_pixels = [cell_bodies(1,:),cell_bodies(end,:), ...
%         cell_bodies(:,1)',cell_bodies(:,end)'];
%     edge_labels = nonzeros(unique(edge_pixels));
%     not_edge_labels = nonzeros(setdiff(unique(cell_bodies(:)),edge_labels));
%     cell_bodies = ismember(cell_bodies,not_edge_labels).*cell_bodies;
    
    %Filter small objects
    cell_props = regionprops(cell_bodies,'Area');
    cell_bodies = ismember(cell_bodies,find([cell_props.Area] > 20)).*cell_bodies;
    
    % Renumber cells to be sequential
    cell_nums = unique(cell_bodies);
    assert(cell_nums(1) == 0, 'Background pixels not found after building object label matrix')
    for i = 2:length(cell_nums)
        cell_bodies(cell_bodies == cell_nums(i)) = i - 1;
    end
    
    % Build cell perimeters image
    cell_bodies_perim = zeros(size(cell_bodies));
    for i = 1:max(cell_bodies(:))
        assert(any(any(cell_bodies == i)), 'Error: can''t find object number %d', i);
        this_ad = zeros(size(cell_bodies));
        this_ad(cell_bodies == i) = 1;
        cell_bodies_perim(bwperim(this_ad)) = i;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Results Output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    output_dir = sprintf('../results/kiehart_testing/02/images/%05d',i_num);
    if (not(exist(output_dir,'dir')))
        mkdir(output_dir);
    end
    
    image_norm = (image - min(image(:)))/range(image(:));
    highlighted_image = create_highlighted_image(image_norm,cell_bodies,'mix_percent',0.5);
    imwrite(highlighted_image,fullfile(output_dir,'highlighted.png'));
    
    imwrite(uint16(cell_bodies),fullfile(output_dir,filenames.objects));
    imwrite(uint16(cell_bodies_perim),fullfile(output_dir,filenames.objects_perim));
    imwrite(uint16(image),fullfile(output_dir,filenames.raw_image));
end

toc;