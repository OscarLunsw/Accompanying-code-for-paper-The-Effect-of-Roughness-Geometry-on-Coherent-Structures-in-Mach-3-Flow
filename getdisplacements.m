%%
% get displacements

function [Dx,Dy] = getdisplacements(dispIm,refIm,WinSize,StepSize,TikonovStrength,LargeDisp,HighRes)
arguments
    dispIm
    refIm
    WinSize double = 32; %window size (px)
    StepSize double = 16; %stepsize (px)
    TikonovStrength double = 1; %fitting coeficant
    LargeDisp logical = false; %large displacement 
    HighRes logical = false; %second pass
end

% deal with empty inputs cause matlab is a lil pissbaby
if isempty(WinSize)
    WinSize = 32;
end
if isempty(StepSize)
    StepSize = 16;
end
if isempty(TikonovStrength)
    TikonovStrength = 1;
end
if isempty(LargeDisp)
    LargeDisp = false;
end
if isempty(HighRes)
    HighRes = false;
end

%% Cross corilation
% compute the size of the discplacement maps based on the window size and the overlap
ImSize=size(refIm);
NumberWindowsH=floor( ( ImSize(2)-(WinSize - StepSize)-4) / StepSize); % keep 4 pixels margin (2 each side) in case of large displacements
NumberWindowsV=floor( ( ImSize(1)-(WinSize - StepSize)-4) / StepSize); % keep 4 pixels margin (2 each side) in case of large displacements

%prealocate variables
Dx=zeros(NumberWindowsV,NumberWindowsH);
Dy=zeros(NumberWindowsV,NumberWindowsH);
DispSize=size(Dx);

Ymax=size(dispIm,1);
Xmax=size(dispIm,2);
jmax=DispSize(1);
imax=DispSize(2);


for i=1:imax
    for j=1:jmax
        % extract the samples of the reference and the experimental image based on the window size
        win1=refIm( (j-1)*StepSize+3:(j-1)*StepSize+WinSize+2 , (i-1)*StepSize+3:(i-1)*StepSize+WinSize+2 );
        win2=dispIm( (j-1)*StepSize+3:(j-1)*StepSize+WinSize+2 , (i-1)*StepSize+3:(i-1)*StepSize+WinSize+2 );
        
        Ind=find(win1==0); % check if the window contains zeros, thus is excluded by the mask
        if length(Ind)>0.2*WinSize*WinSize
            IsExcluded=1;
        else
            IsExcluded=0;
        end
        
        if IsExcluded==0
            if LargeDisp==true
                % estimate interger pixel displacement with correlation
                Crop_size=round(WinSize/4); % evaluate only half the window to save time
                win1_crop=win1(Crop_size:3*Crop_size,Crop_size:3*Crop_size);
                win1_crop=win1_crop-mean(win1_crop(:));
                win2_crop=win2(Crop_size:3*Crop_size,Crop_size:3*Crop_size);
                win2_crop=win2_crop-mean(win2_crop(:));
                WinCorr=conv2(win1_crop,rot90(win2_crop,2),'same');
                
                MaxCorr=max(WinCorr,[],[1 2]);
                [DY,DX]=find(WinCorr==MaxCorr);
                DY=DY(1); % if correlation between two uniform area, DY and DX can have multiple value, so it is needed to chose only one
                DX=DX(1);
                DY=DY-Crop_size-1;
                DX=DX-Crop_size-1;
                % check is the new win2 is conatined in the ExpIm, if not, cancel the displacements
                if  (j-1)*StepSize+3-DY<1 || (i-1)*StepSize+3-DX<1 || (j-1)*StepSize+WinSize+2-DY>Ymax || (i-1)*StepSize+WinSize+2-DX>Xmax
                    DX=0;
                    DY=0;
                end
                % extract a new win2 based on interger pixel shift obtained from correlations
                win2=dispIm( (j-1)*StepSize+3-DY:(j-1)*StepSize+WinSize+2-DY , (i-1)*StepSize+3-DX:(i-1)*StepSize+WinSize+2-DX );
            else
                DX=0;
                DY=0;
            end
            % estimate the displacements from least square method simplified implementation of Appendix A in Pan, Aundi, Xie and Gao, Optics and Laser Engineering, 47 865-874 (2009)
            % define the kernel for computing gradients in the image
            GradVec=[1/12 -2/3 0 2/3 -1/12];
            crop_w=2; % set the number of points to crop at the edge of the subset
  
            win2_x=conv2(win2,GradVec,'same');
            win2_y=conv2(win2,GradVec','same');
            % crop the windows to remove artefact on the edges introduced during gradient computing
            win2_x_crop=win2_x(1+crop_w:end-crop_w,1+crop_w:end-crop_w);
            win2_y_crop=win2_y(1+crop_w:end-crop_w,1+crop_w:end-crop_w);
            win1_crop=win1(1+crop_w:end-crop_w,1+crop_w:end-crop_w);
            win2_crop=win2(1+crop_w:end-crop_w,1+crop_w:end-crop_w);
            
            A=zeros(numel(win2_crop),4);
            for k=1:numel(win2_crop)
                A(k,:)=[ win2_x_crop(k)  win2_y_crop(k)  -win1_crop(k) -1];
            end
            AtA = A'*A;
            lambda = 1e-6 * TikonovStrength * trace(AtA);
            p=(A'*A+ lambda*eye(size(AtA)) )\A'*-win2_crop(:);  % inv(A'*A)*A'*-win2_crop(:);
            dx=p(1);
            dy=p(2);
            if isnan(dx) || isnan(dy)
                dx=0;
                dy=0;
            end
            
            % second pass if enabled
            if HighRes==true
                % interpolate the windows and the gradients so as to take into account first pass displacements
                Win2_Interp=win2;
                dx=-dx;
                dy=-dy;
                if dx>0
                    d_x=dx;
                    ShiftX=0;
                else
                    d_x=1+dx;
                    ShiftX=-1;
                end
                if dy>0
                    d_y=dy;
                    ShiftY=0;
                else
                    d_y=1+dy;
                    ShiftY=-1;
                end
                for m=2:WinSize-1
                    for n=2:WinSize-1
                        pix_center=win2(n+ShiftY,m+ShiftX);
                        pix_N=win2(n+1+ShiftY,m+ShiftX);
                        pix_E=win2(n+ShiftY,m+1+ShiftX);
                        pix_NE=win2(n+1+ShiftY,m+1+ShiftX);
                        Win2_Interp(n,m)=pix_center + d_x*(pix_E-pix_center) + d_y*(pix_N-pix_center) + d_x*d_y*( pix_center+pix_NE - pix_E - pix_N) ;
                    end
                end
                % compute the gradients of the interpolated image
                Win2_x_Interp=conv2(Win2_Interp,GradVec,'same');
                Win2_y_Interp=conv2(Win2_Interp,GradVec','same');
                % crop the windows to remove artefact on the edges introduced during gradient computing
                win2_x_crop=Win2_x_Interp(1+crop_w:end-crop_w,1+crop_w:end-crop_w);
                win2_y_crop=Win2_y_Interp(1+crop_w:end-crop_w,1+crop_w:end-crop_w);
                win2_crop=Win2_Interp(1+crop_w:end-crop_w,1+crop_w:end-crop_w);
                % perform Least Square analysis
                A=zeros(numel(win2_crop),4);
                for k=1:numel(win2_crop)
                    A(k,:)=[ win2_x_crop(k)  win2_y_crop(k)  -win1_crop(k) -1];
                end
                AtA = A'*A;
                lambda = 1e-6 * TikonovStrength * trace(AtA);
                p=(A'*A + lambda*eye(size(AtA)) )\A'*-win2_crop(:);  % inv(A'*A)*A'*-win2_crop(:);
                ddx=p(1);
                ddy=p(2);
                if isnan(ddx) || isnan(ddy)
                    ddx=0;
                    ddy=0;
                end
                dx=-dx+ddx;
                dy=-dy+ddy;
            end
            Dx(j,i)=dx+DX;
            Dy(j,i)=dy+DY;
        else
            Dx(j,i)=0;
            Dy(j,i)=0;
        end
    end
end
end