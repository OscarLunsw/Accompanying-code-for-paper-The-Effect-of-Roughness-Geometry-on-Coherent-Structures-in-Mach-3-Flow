%% Image normalization code
%
function [RefIm,vid] = imagenorm(RefIm1,vid1,ImThresLow,ImThresHigh)
arguments
    RefIm1
    vid1
    ImThresLow double = 0;
    ImThresHigh double = 1;
end

% deal with empty args cause fucking matlab
if isempty(ImThresLow)
    ImThresLow = 0;
end
if isempty(ImThresHigh)
    ImThresHigh = 1;
end

% actually normalize video/image
vid = (vid1-ImThresLow)/(ImThresHigh-ImThresLow);
RefIm = (RefIm1-ImThresLow)/(ImThresHigh-ImThresLow);
end