function setup_experiment(data_folder,varargin)
% SEGMENT_CELLS    Segments cell bodies from E-cadherin images.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
i_p = inputParser;
i_p.StructExpand = true;
i_p.addRequired('data_folder',@(x)exist(x,'dir') == 7);

i_p.parse(data_folder,varargin{:});

filenames = add_filenames_to_struct(struct());

results_folder = regexprep(data_folder,'\/data\/','/results/');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_filenames = dir(data_folder);
image_filenames = image_filenames(3:end);

for i_num = 1:length(image_filenames)
    image = double(imread(fullfile(data_folder,image_filenames(i_num).name)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Results Output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    output_dir = sprintf('%s/images/%05d',results_folder,i_num);
    if (not(exist(output_dir,'dir')))
        mkdir(output_dir);
    end
    
    imwrite(uint16(image),fullfile(output_dir,filenames.raw_image));
    
    image_test = imread(fullfile(output_dir,filenames.raw_image));
    assert(all(all(image == image_test)));
end

toc;