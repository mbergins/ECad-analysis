function track_cells(exp_dir,varargin)

start_time = tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('../find_cell_features'));
filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_dir = fullfile(exp_dir,'images');

image_dirs = dir(base_dir);
image_dirs = image_dirs(3:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image Reading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

objects = cell(0);
tracking_props = cell(0);

image_reading_start = tic;
for i_num=1:length(image_dirs)
    objects{i_num} = imread(fullfile(base_dir,image_dirs(i_num).name,filenames.objects));
        
    for ad_num = 1:max(objects{i_num}(:))
        tracking_props{i_num}(ad_num).assigned = 0;
        tracking_props{i_num}(ad_num).next_obj = [];
    end
end
fprintf('Reading images took ~%d sminutes.\n', round(toc(image_reading_start)/60));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Object Assocation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assign_start = tic;
for i_num=1:(length(objects) - 1)
    pix_sim = calc_pix_sim(objects{i_num},objects{i_num+1});
    cent_dist = calc_cent_dist(objects{i_num},objects{i_num+1});
        
    %Start by searching for reciprical high pixel similarity matches,
    %defined as those objects that overlap a single object in the next frame by
    %50% or more
    [start_obj_hits, end_obj_hits] = find(pix_sim > 0.5);    
    for i = 1:length(start_obj_hits)
        start_obj = start_obj_hits(i);
        end_obj = end_obj_hits(i);
        
        %check to make sure these two objects are unique in their lists, if
        %so, make the connection and block the row (start_obj) and column
        %(end_obj)
        if (sum(start_obj == start_obj_hits) == 1 && ...
            sum(end_obj == end_obj_hits) == 1)
            
            tracking_props{i_num}(start_obj).next_obj = end_obj;
            
            pix_sim(start_obj,:) = NaN;
            pix_sim(:,end_obj) = NaN;

            cent_dist(start_obj,:) = NaN;
            cent_dist(:,end_obj) = NaN;
        end
    end
    
    %This loop finds any remaining pixel similarity measures above 20% and
    %matches them, one-by-one to their corresponding end objects. This code
    %mostly helps clear out the objects that overlap two objects with fairly
    %high similarity, but not one single object completely.
    while (any(any(pix_sim > 0.2)))
        [start_obj,end_obj] = find(pix_sim == max(pix_sim(:)),1,'first');
        
        tracking_props{i_num}(start_obj).next_obj = end_obj;
        
        pix_sim(start_obj,:) = NaN;
        pix_sim(:,end_obj) = NaN;

        cent_dist(start_obj,:) = NaN;
        cent_dist(:,end_obj) = NaN;
    end

    %This section of the code would use the distance between the centroids
    %to fill out the rest of the tracking list, but it isn't needed for the
    %puncta tracking
%     [start_obj_hits, end_obj_hits] = find(cent_dist < 10);
%     for i = 1:length(start_obj_hits)
%         start_obj = start_obj_hits(i);
%         end_obj = end_obj_hits(i);
%         
%         %check to make sure these two objects are unique in their lists, if
%         %so, make the connection and block the row (start object) and column
%         %(end object)
%         if (sum(start_obj == start_obj_hits) == 1 && ...
%             sum(end_obj == end_obj_hits) == 1)
%             
%             tracking_props{i_num}(start_obj).next_obj = end_obj;
%             
%             pix_sim(start_obj,:) = NaN;
%             pix_sim(:,end_obj) = NaN;
% 
%             cent_dist(start_obj,:) = NaN;
%             cent_dist(:,end_obj) = NaN;
%         end
%     end
    
    if (mod(i_num,10)==0)
        runtime_now = toc(assign_start);
        estimated_remaining = round(((runtime_now/i_num)*(length(objects) - 1 - i_num))/60);
        fprintf('Done with image %d, estimating %d minutes left.\n',i_num,estimated_remaining);
    end
end

tracking_build_start = tic;
tracking_mat = convert_tracking_props_to_matrix(tracking_props);

%Dealing with the case where the last frames of a movie don't contain any
%objects, so the tracking matrix is truncated, resulting in shortened
%visualizations. In this case, pad the tracking mat with -1 values to
%ensure that all frames are represented, even if empty.
if (size(tracking_mat,2) < length(image_dirs))
    padding_mat = -1 * ones(size(tracking_mat,1),length(image_dirs) - size(tracking_mat,2));
    tracking_mat = [tracking_mat,padding_mat];
end

tracking_build_time = round(toc(tracking_build_start)/60);
fprintf('Tracking matrix building took %d minutes.\n',tracking_build_time);

if (not(exist(fullfile(base_dir,'..','tracking_matrices'),'dir')))
    mkdir(fullfile(base_dir,'..','tracking_matrices'))
end
toc(start_time);

output_file = fullfile(base_dir, image_dirs(1).name,filenames.tracking);

%If the tracking matrix is empty, don't output anything and don't make a
%folder for the empty matrix
if (not(any(size(tracking_mat) == 0)))
    %check each column of the tracking matrix to make sure each object number
    %occurs once and only once per column
    for col_num = 1:size(tracking_mat,2)
        all_tracking_nums = sort(tracking_mat(:,col_num));
        for object_num = 1:max(all_tracking_nums)
            assert(sum(all_tracking_nums == object_num) == 1);
        end
    end
    
    if (not(exist(fileparts(output_file),'dir'))), mkdir(fileparts(output_file)); end
    
    csvwrite(fullfile(base_dir, image_dirs(1).name,filenames.tracking),tracking_mat);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similarity Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pix_sim = calc_pix_sim(ads_1,ads_2)
    pix_sim = zeros(max(ads_1(:)), max(ads_2(:)));
    
    for start_ad=1:max(ads_1(:))
        this_ad = ads_1 == start_ad;
        
        ad_size = sum(sum(this_ad));
        
        overlap_pix = ads_2(this_ad);
        
        overlap_pix = overlap_pix(overlap_pix ~= 0);
        if (isempty(overlap_pix))
            continue;
        end
        
        unique_overlap = unique(overlap_pix);
        
        for end_ad = unique_overlap'
            pix_sim(start_ad,end_ad) = sum(overlap_pix == end_ad)/ad_size;
        end
    end
end

function cent_dist = calc_cent_dist(ads_1,ads_2)
    props_1 = regionprops(ads_1,'Centroid');
    props_2 = regionprops(ads_2,'Centroid');
    
    centroid_1 = reshape([props_1.Centroid],2,[])';
    centroid_2 = reshape([props_2.Centroid],2,[])';
    cent_dist = pdist2(centroid_1,centroid_2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tracking Matrix Production
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tracking_matrix = convert_tracking_props_to_matrix(tracking_props)

objects = struct('start',{},'sequence',{});
tracking_num = 1;

[i_num,obj_num] = find_unassigned_obj(tracking_props);

while (i_num ~= 0 && obj_num ~= 0)
    if (length(objects) < tracking_num)
        objects(tracking_num).start = i_num;
    end
    
    objects(tracking_num).sequence = [objects(tracking_num).sequence, obj_num];
    tracking_props{i_num}(obj_num).assigned = 1;
    
    %pick out the next object to follow
    [i_num, obj_num] = follow_to_next_obj(tracking_props,i_num,obj_num);
    
    if (i_num == 0)
        [i_num, obj_num] = find_unassigned_obj(tracking_props);
        tracking_num = tracking_num + 1;
        if (mod(tracking_num,1000) == 0)
            fprintf('Done with building tracking for %d objects.\n', tracking_num);
        end
    end
end

tracking_matrix = zeros(length(objects),length(tracking_props));

for obj_num = 1:length(objects)
    col_range = objects(obj_num).start:(objects(obj_num).start + length(objects(obj_num).sequence) - 1);
    tracking_matrix(obj_num,col_range) = objects(obj_num).sequence;
end

for col_num = 1:size(tracking_matrix,2)
    this_col = tracking_matrix(:,col_num);
    this_col = this_col(this_col ~= 0);
    
    this_col = sort(unique(this_col))';
    
    %empty columns mean there weren't any objects in that time step, so
    %check for that in the following assert first, then if there were
    %objects, make sure all were accounted for
    assert(isempty(this_col) || all(this_col == 1:max(this_col)))
    
    assert((isempty(tracking_props{col_num}) && isempty(this_col)) || ...
        (length(tracking_props{col_num}) == max(this_col)));
end

end

function [i_num,obj_num] = find_unassigned_obj(tracking_props)

%scan through the tracking props data to find an entry where the assigned
%value is false, then immediatly return, taking the current value of i_num
%and obj_num with the return
for i_num=1:length(tracking_props)
    for obj_num=1:max(size(tracking_props{i_num}))
        try
            if (tracking_props{i_num}(obj_num).assigned)
                continue;
            else
                return;
            end
        catch
            continue;
        end
    end
    tracking_props{i_num} = [];
end

%we won't get to this part of the code unless there aren't unassigned
%objects left, in that case return the exit signal
i_num = 0;
obj_num = 0;

return;

end

function [i_num,obj_num] = follow_to_next_obj(tracking_props,i_num,obj_num)

try %#ok<*TRYNC>
    if (size(tracking_props{i_num}(obj_num).next_obj,2) == 0)
        i_num = 0;
        obj_num = 0;
        return;
    else
        obj_num = tracking_props{i_num}(obj_num).next_obj;
        i_num = i_num + 1;
        return;
    end
end

i_num = 0;
obj_num = 0;

return;

end
