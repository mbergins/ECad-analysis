function data_set = read_in_file_set(this_dir,filenames)
data_set = struct;

data_set.this_dir = this_dir;

data_set.image  = double(imread(fullfile(this_dir, filenames.raw_image)));
data_set.image_norm = (data_set.image - min(data_set.image(:)))/range(data_set.image(:));

data_set.labeled_cells = imread(fullfile(this_dir, filenames.objects));
data_set.labeled_cells_perim = imread(fullfile(this_dir, filenames.objects_perim));