%%
% load image
function RefIm = loadimage(filename)
    RefIm = im2double(imread(filename));
    %converts to greyscale if image is RGB
    if size(RefIm,3) == 3
        RefIm = rgb2gray(RefIm);
    end
disp("Image Loaded")
end