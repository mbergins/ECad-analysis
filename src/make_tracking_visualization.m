function make_tracking_visualization(exp_dir,varargin)

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

image_dir = fullfile(exp_dir,'images');
single_image_dirs = dir(image_dir);

single_image_dirs = single_image_dirs(3:end);

tracking_file = fullfile(image_dir, single_image_dirs(1).name,filenames.tracking);
if (not(exist(tracking_file,'file')))
    tracking_matrix = zeros(1,size(single_image_dirs,1));
else
    tracking_matrix = csvread(tracking_file);
end

lineage_cmap = jet(size(tracking_matrix,1));
% lineage_cmap = lineage_cmap(randperm(size(tracking_matrix,1)),:);

lineage_cmap = assign_lineage_to_cmap(tracking_matrix);

convert_avail = not(system('which convert > /dev/null'));
for i_num = 1:length(single_image_dirs)
    current_data = read_in_file_set(fullfile(image_dir,single_image_dirs(i_num).name),filenames);
    
    %thicken the cell perimeter, easier for visualization
%     labeled_cells_perim_thick = thicken_perimeter(current_data.labeled_cells_perim,current_data.labeled_cells);
    labeled_cells_perim_thick = current_data.labeled_cells_perim;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Image Creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    live_rows = tracking_matrix(:,i_num) > 0;
    ad_nums = tracking_matrix(live_rows,i_num);
    tracking_nums = find(live_rows);
    
    %Build the unique lineage highlighted image
    this_cmap(ad_nums,:) = lineage_cmap(live_rows,:); %#ok<AGROW>
    highlighted_image = create_highlighted_image(current_data.image_norm, ...
        labeled_cells_perim_thick,'color_map',this_cmap);
    
    output_folder = fullfile(exp_dir,'visualizations','tracking');
    if (not(exist(output_folder,'dir')))
        mkdir(output_folder);
    end
    output_file = fullfile(output_folder,sprintf('%05d.png',i_num));
    
    imwrite(highlighted_image,output_file);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Image Labeling
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (not(convert_avail) || not(exist(tracking_file,'file'))), continue; end
    props = regionprops(current_data.labeled_cells,'Centroid');
    
    i_size = size(highlighted_image);
    for cell_num=1:length(props)
        if (props(cell_num).Centroid(1) > i_size(2)*0.95)
            props(cell_num).Centroid(1) = i_size(2)*0.95;
        end
        if (props(cell_num).Centroid(2) > i_size(1)*0.95)
            props(cell_num).Centroid(2) = i_size(1)*0.95;
        end
    end
    
    all_annotate = '';
    for cell_num = 1:length(ad_nums)
        pos_str = ['+',num2str(props(ad_nums(cell_num)).Centroid(1)),'+',num2str(props(ad_nums(cell_num)).Centroid(2))];
        label_str = [' "',num2str(tracking_nums(cell_num)),'"'];
        all_annotate = [all_annotate, ' -annotate ', pos_str, label_str, ' ']; %#ok<AGROW>
    end
    
    rgb_code = '';
    if (mean(current_data.image_norm(:)) > 0.5)
        rgb_code = '''rgb(0,0,0)'' ';        
    else
        rgb_code = '''rgb(255,255,255)'' ';
    end
    command_str = ['convert ', output_file, ' -font VeraBd.ttf -pointsize 10 -fill ', ...
        rgb_code,all_annotate, output_file, '; '];
    system(command_str);
end

disp(['Done with ', exp_dir]);
toc;
end


function lineage_cmap = assign_lineage_to_cmap(tracking_matrix)

cell_counts = tracking_matrix > 0;
max_cell_number = max(sum(tracking_matrix > 0));

cmap = jet(max_cell_number);
lineage_to_cmap = zeros(size(tracking_matrix,1),1);

for i_num = 1:size(tracking_matrix,2)
    for j = 1:size(tracking_matrix,1)
        if (tracking_matrix(j,i_num) <= 0), continue; end

        %Unique lineage colors
        if (lineage_to_cmap(j) == 0)
            used_c_nums = sort(lineage_to_cmap(tracking_matrix(:,i_num) > 0));
            used_c_nums = used_c_nums(used_c_nums ~= 0);

            taken_nums = zeros(1,max_cell_number);
            taken_nums(used_c_nums) = 1;
            taken_dists = bwdist(taken_nums);

            try
                lineage_to_cmap(j) = find(taken_dists == max(taken_dists),1,'first');
            catch map_error %#ok<NASGU>
                assert(~any(taken_dists == max(taken_dists)), 'Error: could not find a possible color number in image number %d',padded_i_num);
            end
        end
    end
end

lineage_cmap = cmap(lineage_to_cmap,:);

end