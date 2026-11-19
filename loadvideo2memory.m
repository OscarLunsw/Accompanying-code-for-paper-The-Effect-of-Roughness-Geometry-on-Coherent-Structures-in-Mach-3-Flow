%%
% loads video to an array
% If Fstart and Fend are specified the frames between Fstart and Fend
% inclusive will be loaded otherwise the whole video will be loaded

function [N,vid]=loadvideo2memory(filename,Fstart,Fend)
arguments
    filename string
    Fstart double = 1
    Fend double = inf
end
% deal with empty arguments cause matlab is a little baby
if isempty(Fstart)
    Fstart = 1;
end
if isempty(Fend)
    Fend = inf;
end

%load video
vid_file = VideoReader(filename);

% find number of frames in video
if isinf(Fend)
    Fend = floor(vid_file.Duration * vid_file.FrameRate);
end
N = Fend - Fstart + 1;

%prealocate vid for speed
vid = zeros(vid_file.Height,vid_file.Width,N,'double');

%loop through frames extracting them to vid
i=1;
while Fstart+i-1 <= Fend
    fprintf('Loading Frame %g \n' ,i) %print status to workspace
    frame = im2double(read(vid_file,Fstart+i-1));
    
    % converts to greyscale if video is RGB
    if size(frame,3) == 3
        frame = rgb2gray(frame);
    end
    vid(:,:,i) = frame;
    i=i+1;
end
disp("Video Load Finished")