% writes the displacements to an .avi file

function writeavi(mag,N,outfilename,framerate,lower_cutoff,upper_cutoff)
    arguments
        mag double
        N int64
        outfilename string
        framerate double = 1
        lower_cutoff double = []
        upper_cutoff double = []
    end
    % Default cutoffs to mat2gray behavior
    if isempty(lower_cutoff)
        lower_cutoff = min(mag);
    end
    if isempty(upper_cutoff)
        upper_cutoff = max(mag);
    end

    % Setup output file
    outVid = VideoWriter(outfilename);
    outVid.FrameRate = framerate; 
    outVid.Quality = 100;
    open(outVid);
    disp("Writing Video")

    % Write output file
    mag_01 = mat2gray(mag,[lower_cutoff,upper_cutoff]);

    for i = 1:N
        frame = im2uint8(mag_01(:,:,i));
        writeVideo(outVid,frame);
    end

    close(outVid);
    disp("Done!")
end