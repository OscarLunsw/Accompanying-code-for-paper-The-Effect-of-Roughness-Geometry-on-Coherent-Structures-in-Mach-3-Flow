%% Background Oriented Schlieren (BOS) Processing Script

% Registers distorted and reference BOS images
% Calculates displacement fields
% Removes outliers
% Reconstructs refractive index and density fields
% Estimates boundary layer velocity profile using Crocco-Busemann relation

% Author: Oscar Lewis
% Date: 06/08/2026

clear; clc; close all;

%% INPUT PARAMETERS

% File paths
% Add required function folders
addpath(genpath('C:\Users\oscar\Documents\UNSW\Thesis\BOS'))

% Working directory containing images and output files
input.folder = fullfile( ...
    'C:\Users\oscar\Documents\UNSW\Thesis\BOS', ...
    'Final testing 23-07', ...
    'Tirangles');

% BOS images
input.refFile  = fullfile(input.folder,'Ref1.tif');
input.dispFile = fullfile(input.folder,'T625.tif');

% Image geometry
% Image height corresponding to tunnel insert/floor location (pixels)
input.floor_pix = 4646;

% Camera pixel dimensions
cal.pixSize = 2.5e-6;       % Camera pixel size [m]

% Total image width
input.nCols = 5120;

% BOS calibration parameters
cal.G       = 2.26261e-4;   % Gladstone-Dale constant [m^3/kg]
cal.Z_T     = 0.101;        % Optical path length / tunnel width [m]
cal.factor  = 0.187;        % Z_i/(Z_A + Z_T + Z_window)

% BOS processing parameters
% Displacement calculation
process.WinSize          = 32;      % Correlation interrogation size [px]
process.StepSize         = 16;      % Correlation spacing [px]
process.TikhonovStrength = 1;       % Gradient reconstruction regularisation
process.LargeDisp        = true;    % Enable large displacement mode
process.HighRes          = true;    % Enable second high-resolution pass

% Outlier filtering
process.outlierThreshold = 1;       % Maximum displacement deviation [px]

% Velocity reconstruction parameters
vel.R  = 287;      % Specific gas constant [J/kg/K]
vel.cp = 1005;     % Specific heat capacity [J/kg/K]

%% Tunnel flow conditions
load(fullfile(input.folder,"DAQ_data/","T625.mat"))

flow.P_s      = S.testConditions.trigger.Ps_kPa;      % Static pressure [kPa]
flow.P_t      = S.testConditions.trigger.Pt_kPa;     % Stagnation pressure [kPa]

flow.T_inf    = S.testConditions.trigger.Tinf_K;     % Static temperature [K]
flow.T_0      = S.testConditions.trigger.T0_K;     % Stagnation temperature [K]

flow.Mach     = S.testConditions.trigger.Mach;
flow.U_inf    = S.testConditions.trigger.Uinf_mps;     % Freestream velocity [m/s]
flow.rho_inf  = S.testConditions.trigger.rhoinf_kgm3;    % Freestream density [kg/m^3]
clear S

load(fullfile(input.folder,"DAQ_data/","Pressure_tap_data.mat"))
flow.P_surf   = S.testConditions.trigger.Psurf_kpa;     % Surface Pressure [kPa]
clear S

%% PREPROCESSING

disp('Preparing image geometry')

% Freestream refractive index
cal.n0 = cal.G*flow.rho_inf + 1;

% Determine boundary region used for velocity extraction
input.boundary_pix = floor( ...
    input.floor_pix - ...
    0.012*cal.factor/cal.pixSize);

% Pixel coordinate vector
x = 1:input.nCols;

% Define masking boundaries
% Region between these curves contains unwanted tunnel geometry
yTop = 0.02*input.floor_pix + 0.45*x;
yBot = 0.50*input.floor_pix + 0.30*x;

% Generate polygon mask
xPoly = [x fliplr(x)];
yPoly = [yTop fliplr(yBot)];

mask = ~poly2mask( ...
    xPoly, ...
    yPoly, ...
    input.floor_pix, ...
    input.nCols);

% Remove region below boundary layer
mask(input.boundary_pix:end,:) = false;


%% LOAD IMAGES

disp('Loading images')

% Crop images to remove tunnel insert region
refIm = loadimage(input.refFile);
refIm = refIm(1:input.floor_pix,:);

dispIm = loadimage(input.dispFile);
dispIm = dispIm(1:input.floor_pix,:);

%% IMAGE REGISTRATION

disp('Aligning images')

% Apply mask during registration only
refMasked  = refIm .* mask;
dispMasked = dispIm .* mask;

% Configure registration optimizer and similarity metric
[registration.optimizer, registration.metric] = ...
    imregconfig('monomodal');

% Adjust optimizer settings if required
registration.optimizer.MaximumIterations = 300;

% Determine image transformation
registration.tform = imregtform( ...
    dispMasked, ...
    refMasked, ...
    'translation', ...
    registration.optimizer, ...
    registration.metric);

% Apply transformation to original distorted image
dispIm = imwarp( ...
    dispIm, ...
    registration.tform, ...
    'OutputView', ...
    imref2d(size(refIm)));

%% DISPLACEMENT CALCULATION

disp('Calculating displacement fields')

[Dx_raw,Dy_raw] = getdisplacements( ...
    dispIm, ...
    refIm, ...
    process.WinSize, ...
    process.StepSize, ...
    process.TikhonovStrength, ...
    process.LargeDisp, ...
    process.HighRes);

save(fullfile(input.folder,'Dx_Dy_raw.mat'), ...
    'Dx_raw','Dy_raw');

clear dispIm refIm

%% OUTLIER REMOVAL

disp('Removing displacement outliers')

% Local median displacement
Dx_med = medfilt2(Dx_raw,[3 3]);
Dy_med = medfilt2(Dy_raw,[3 3]);

% Identify outliers
Dx_outliers = abs(Dx_raw-Dx_med) > process.outlierThreshold;
Dy_outliers = abs(Dy_raw-Dy_med) > process.outlierThreshold;

% Replace with local median
Dx_filtered = Dx_raw;
Dy_filtered = Dy_raw;

Dx_filtered(Dx_outliers) = Dx_med(Dx_outliers);
Dy_filtered(Dy_outliers) = Dy_med(Dy_outliers);

save(fullfile(input.folder,'Dx_Dy_filtered.mat'), ...
    'Dx_filtered','Dy_filtered');

%% DENSITY RECONSTRUCTION
% Gladstone-Dale Relation

disp('Calculating density field')

% Convert displacement into angular deflection
cal.theta_x = ...
    (Dx_filtered .* cal.pixSize) ./ ...
    (cal.factor*cal.Z_T);

cal.theta_y = ...
    (Dy_filtered .* cal.pixSize) ./ ...
    (cal.factor*cal.Z_T);

% Calculate refractive index gradients
cal.n_x = cal.theta_x/cal.Z_T;
cal.n_y = cal.theta_y/cal.Z_T;

% Integrate gradients to obtain refractive index field
pixelSpacing = ...
    cal.pixSize*process.StepSize/cal.factor;

cal.n = intgrad2( ...
    cal.n_x,...
    cal.n_y,...
    pixelSpacing,...
    pixelSpacing);

% Correct arbitrary integration offset using freestream region
cal.n_offset = ...
    cal.n0 - mean(cal.n(1:5,end-5:end),"all");

cal.n = cal.n + cal.n_offset;

% Convert refractive index to density
cal.rho = (cal.n-1)/cal.G;

save(fullfile(input.folder,'Density.mat'), ...
    '-struct','cal');

%% VELOCITY PROFILE RECONSTRUCTION
% Crocco-Busemann Relation

disp('Calculating velocity profile')

% Convert floor location to displacement grid coordinates
vel.floor_pix = ...
    floor((input.floor_pix-process.WinSize)/process.StepSize)+1;

vel.boundary_pix = ...
    floor((input.boundary_pix-process.WinSize)/process.StepSize)+1;

% Extract boundary layer density
vel.rho = ...
    cal.rho(vel.boundary_pix:vel.floor_pix,:);

% Spatial average
vel.rho_mean = mean(vel.rho,2);

% Temperature from ideal gas law
vel.temp = ...
    (flow.P_surf*1e3)./(vel.R*vel.rho_mean);

% Crocco-Busemann velocity relation
vel.u = flip( ...
    sqrt(2*vel.cp*(vel.temp(end)-vel.temp)) );

% Wall-normal coordinate
vel.y = ...
    (0:length(vel.u)-1)' * pixelSpacing;

save(fullfile(input.folder,'Boundary_profile.mat'), ...
    '-struct','vel');

figure
plot(vel.y,vel.u)
xlabel('y [m]')
ylabel('Velocity [m/s]')
grid on

savefig(fullfile(input.folder,'boundary_profile'))
close

%% POST PROCESSING

disp('Saving displacement magnitude')

% Displacement magnitude
displacementMag = hypot(Dx_filtered,Dy_filtered);

imwrite( ...
    mat2gray(displacementMag), ...
    fullfile(input.folder,'Pix_disp.png'));

disp('BOS processing complete')
