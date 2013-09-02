function cfg_struct = add_filenames_to_struct(cfg_struct)
%ADD_FILENAMES_TO_STRUCT    Takes in a struct, to which the filenames for
%                          various files made in the pipeline are added

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg_struct.raw_image = 'raw_image.png';

cfg_struct.objects = 'cell_labeled.png';
cfg_struct.objects_perim = 'cell_labeled_perim.png';

cfg_struct.tracking = '../../tracking_matrices/tracking_seq.csv';