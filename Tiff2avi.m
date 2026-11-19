% Folder containing TIFF files
imageFolder = 'C:\Users\oscar\Documents\UNSW\Thesis\BOS\Flame test 2 data\candle_tiff\Candle';

% Get list of TIFF files
imageFiles = dir(fullfile(imageFolder, '*.tif'));

% Sort files alphabetically
[~, idx] = sort({imageFiles.name});
imageFiles = imageFiles(idx);

% Create video writer
outputVideo = VideoWriter('output.avi');
outputVideo.FrameRate = 30;   % Frames per second
open(outputVideo);

% Loop through images
for k = 1:length(imageFiles)

    % Read image
    filename = fullfile(imageFolder, imageFiles(k).name);
    img = imread(filename);

    % Convert grayscale to RGB (optional but recommended for AVI compatibility)
    if size(img,3) == 3
        img = rgb2gray(img);
    end
    
    img = im2uint8(mat2gray(img));
   
    % Write frame
    writeVideo(outputVideo, img);

end

% Close video file
close(outputVideo);

disp('Video creation complete.');