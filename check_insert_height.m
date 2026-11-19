%% check insert height

refIm = loadimage('Ref1.tif');
refIm = mat2gray(refIm);
h = imshow(refIm);
hp = impixelinfo;
set(hp,"Position",[5 1 300 20]);