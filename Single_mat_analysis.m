% This file is a companion script to CONVERTER_bin_to_mat.
%
% It loads a converted wind tunnel DAQ .mat file, examines tunnel
% conditions over time, detects the camera trigger, calculates Mach number,
% and stores the test conditions at the trigger time back into the struct.
%
% Input:
%   .mat file containing struct S from DAQ converter
%
% Output:
%   - pressure, temperature, trigger-time and Mach plots
%   - updated S.analysis and S.testConditions fields
%
% Written by Erik Mueller, with help from ChatGPT.

clear;
clc;

%% USER SETTINGS

[FileName,PathName] = uigetfile();
matFile = fullfile(PathName,FileName);

% Channel IDs
CH_STATIC_PRESSURE = 80;     % Baratron static pressure
CH_TOTAL_PRESSURE  = 81;     % Wika stagnation/total pressure
CH_SURFACE_PRESSURE = 83;    % Surface Pressure Tap
CH_TRIGGER         = 113;    % camera trigger witness

% Gas constants
gamma = 1.4;
R_air = 287;                 % J/(kg.K)

% Trigger detection settings
triggerThreshold_V = 4.5;    % trigger expected to jump to approximately 5 V
minHighDuration_s  = 1e-3;   % require signal to remain high for 1 ms

% Plot settings
Mach_ylim = [0 3.5];

% Save updated S back to the same .mat file?
saveUpdatedStruct = true;

%% LOAD CONVERTED DAQ FILE

loadedData = load(matFile);

if ~isfield(loadedData, 'S')
    error('Selected .mat file does not contain expected struct S.');
end

S = loadedData.S;
clear loadedData;

printFileSummary(S);

%% EXTRACT COMMON VARIABLES

t = S.time.highRate_s;

channelIDs = [S.channels.id];

Ps_kPa    = getChannelProcessedValue(S, CH_STATIC_PRESSURE, 'static pressure');
Pt_kPa    = getChannelProcessedValue(S, CH_TOTAL_PRESSURE,  'stagnation pressure');
Psurf_kpa = getChannelProcessedValue(S, CH_SURFACE_PRESSURE,  'surface pressure');
trigger_V = getChannelProcessedValue(S, CH_TRIGGER,         'trigger witness');

t_T = S.temperature.time_s;
T0_K = S.temperature.K;       % stagnation temperature, low-rate

%% DETECT TRIGGER TIME

[triggerTime_s, triggerIdx] = detectTriggerTime( ...
    t, trigger_V, triggerThreshold_V, minHighDuration_s, S.acquisition.sampleRate_Hz);

% If no trigger is detected, assume trigger occurred at DAQ start.
% This covers cases where the trigger was manually set high before or at
% the moment acquisition began.
if isempty(triggerTime_s)

    triggerTime_s = 0.000;
    triggerIdx    = 1;

    warning('No valid trigger event detected. Setting trigger time to DAQ start: %.3f s.', ...
        triggerTime_s);

end

fprintf('\nTrigger detection:\n');
fprintf('  Trigger time:          %.9f s\n', triggerTime_s);
fprintf('  Trigger sample index:  %.0f\n', triggerIdx);
fprintf('  Threshold used:        %.3f V\n', triggerThreshold_V);
fprintf('  Minimum high duration: %.6f s\n', minHighDuration_s);

%% CALCULATE MACH NUMBER AND FLOW PROPERTIES

Mach = [];
PressureRatio = [];
T0_high_K = [];
Tinf_K = [];
a_inf_mps = [];
Uinf_mps = [];
rhoinf_kgm3 = [];

if ~isempty(Ps_kPa) && ~isempty(Pt_kPa)

    PressureRatio = Pt_kPa ./ Ps_kPa;

    validPR = PressureRatio >= 1 & isfinite(PressureRatio);

    Mach = NaN(size(PressureRatio));
    Mach(validPR) = sqrt( ...
        (2/(gamma - 1)) .* ...
        (PressureRatio(validPR).^((gamma - 1)/gamma) - 1) );

    if ~isempty(T0_K) && ~isempty(t_T)

        % Interpolate low-rate stagnation temperature onto high-rate time
        T0_high_K = interp1(t_T, T0_K, t, 'linear', 'extrap');

        % Static/freestream temperature
        Tinf_K = T0_high_K ./ (1 + ((gamma - 1)/2).*Mach.^2);

        % Speed of sound and velocity
        a_inf_mps = sqrt(gamma .* R_air .* Tinf_K);
        Uinf_mps  = Mach .* a_inf_mps;

        % Freestream density
        rhoinf_kgm3 = (Ps_kPa .* 1000) ./ (R_air .* Tinf_K);

    else
        warning('Temperature data not found. Mach calculated, but T_inf, U_inf and rho_inf were not calculated.');
    end

else
    warning('Static and/or stagnation pressure data missing. Mach number was not calculated.');
end

%% STORE MACH TIME SERIES IN STRUCT

S.analysis = struct();

S.analysis.Mach.gamma = gamma;
S.analysis.Mach.R_air = R_air;
S.analysis.Mach.time_s = t;
S.analysis.Mach.pressureRatio_Pt_over_Ps = PressureRatio;
S.analysis.Mach.M = Mach;
S.analysis.Mach.assumption = ...
    'Mach number calculated from isentropic Pt/Ps relation using absolute pressures.';

if ~isempty(T0_high_K)
    S.analysis.Mach.T0_highRate_K = T0_high_K;
    S.analysis.Mach.Tinf_K = Tinf_K;
    S.analysis.Mach.a_inf_mps = a_inf_mps;
    S.analysis.Mach.Uinf_mps = Uinf_mps;
    S.analysis.Mach.rhoinf_kgm3 = rhoinf_kgm3;
end

if ~isempty(triggerTime_s)
    S.analysis.trigger.time_s = triggerTime_s;
    S.analysis.trigger.sampleIndex_highRate = triggerIdx;
    S.analysis.trigger.threshold_V = triggerThreshold_V;
    S.analysis.trigger.minHighDuration_s = minHighDuration_s;
end

%% SAVE TEST CONDITIONS AT TRIGGER TIME INTO STRUCT

if ~isempty(triggerTime_s)

    idx = triggerIdx;

    triggerConditions = struct();

    triggerConditions.time_s = triggerTime_s;
    triggerConditions.sampleIndex_highRate = idx;
    triggerConditions.sampleIndex_temperature = nearestIndex(t_T, triggerTime_s);

    triggerConditions.Ps_kPa = valueAt(Ps_kPa, idx);
    triggerConditions.Pt_kPa = valueAt(Pt_kPa, idx);
    triggerConditions.Psurf_kpa = valueAt(Psurf_kpa,idx);

    triggerConditions.Tinf_K = valueAt(Tinf_K, idx);
    triggerConditions.T0_K   = valueAt(T0_high_K, idx);

    triggerConditions.Mach = valueAt(Mach, idx);
    triggerConditions.Uinf_mps = valueAt(Uinf_mps, idx);
    triggerConditions.rhoinf_kgm3 = valueAt(rhoinf_kgm3, idx);

    triggerConditions.gamma = gamma;
    triggerConditions.R_air = R_air;

    triggerConditions.notes = ...
        'Conditions evaluated at nearest high-rate sample to detected trigger time.';

    S.testConditions.trigger = triggerConditions;

    printTriggerConditions(triggerConditions);

else
    warning('Trigger time not available. Test conditions at trigger were not stored.');
end

%% PLOT WIND TUNNEL CONDITIONS AND MACH NUMBER

figure;

% -------------------------------------------------------------------------
% Subplot 1: pressures and stagnation temperature
% -------------------------------------------------------------------------

subplot(2,1,1)
hold on

yyaxis left
plotIfAvailable(t, Ps_kPa, 'P_s', 1.2, "r-");
plotIfAvailable(t, Pt_kPa, 'P_t', 1.2, "g-");
plotIfAvailable(t, Psurf_kpa, 'P_{surf}', 1.2, "b-");
ylabel('Pressure [kPa]')

yyaxis right
plotIfAvailable(t_T, T0_K, 'T_0', 1.2, 'o-');
ylabel('Temperature [K]')

addTriggerLine(triggerTime_s);

grid on
xlim([t(1), t(end)])
legend('Location', 'best')

title(sprintf('Wind Tunnel Conditions: %s', S.file.sourceName), ...
    'Interpreter', 'none')

% -------------------------------------------------------------------------
% Subplot 2: Mach number
% -------------------------------------------------------------------------

subplot(2,1,2)
hold on

plotIfAvailable(t, Mach, 'Mach', 1.2);

addTriggerLine(triggerTime_s);

xlabel('Time [s]')
ylabel('Mach number, M [-]')
grid on
xlim([t(1), t(end)])
ylim(Mach_ylim)
legend('Location', 'best')
title('Mach Number')

%% SAVE UPDATED STRUCT BACK TO MAT FILE

if saveUpdatedStruct
    save(matFile, 'S', '-v7.3');

    fprintf('\nUpdated .mat file saved with analysis and test conditions:\n');
    fprintf('  %s\n', matFile);
end

%% LOCAL FUNCTIONS

function printFileSummary(S)

    fprintf('\nLoaded DAQ file:\n');
    fprintf('  Source file:      %s\n', S.file.sourceFilename);
    fprintf('  Samples:          %.0f\n', S.acquisition.numberOfSamples);
    fprintf('  Sample rate:      %.0f Hz\n', S.acquisition.sampleRate_Hz);
    fprintf('  Run time:         %.3f s\n', S.acquisition.runTime_s);
    fprintf('  Number channels:  %.0f\n', S.acquisition.numberOfChannels);
    fprintf('  Channel IDs:      ');
    fprintf('%.0f ', S.acquisition.channelIDs);
    fprintf('\n\n');

end

function signal = getChannelProcessedValue(S, channelID, signalName)

    channelIDs = [S.channels.id];
    idx = find(channelIDs == channelID, 1);

    if isempty(idx)
        signal = [];
        warning('Channel %.0f for %s not found.', channelID, signalName);
    else
        signal = S.channels(idx).processedValue;
    end

end

function [triggerTime_s, triggerIdx] = detectTriggerTime(t, trigger_V, threshold_V, minHighDuration_s, Fs)

    triggerTime_s = [];
    triggerIdx = [];

    if isempty(trigger_V)
        warning('Trigger signal is empty. Cannot detect trigger time.');
        return
    end

    minHighSamples = max(1, round(minHighDuration_s * Fs));

    aboveThreshold = trigger_V >= threshold_V;
    risingIdx = find(diff([false; aboveThreshold]) == 1);

    for ii = 1:numel(risingIdx)

        idx1 = risingIdx(ii);
        idx2 = min(idx1 + minHighSamples - 1, numel(aboveThreshold));

        if all(aboveThreshold(idx1:idx2))
            triggerIdx = idx1;
            triggerTime_s = t(idx1);
            return
        end

    end

end

function idx = nearestIndex(x, x0)

    if isempty(x) || isempty(x0)
        idx = [];
    else
        [~, idx] = min(abs(x - x0));
    end

end

function y = valueAt(x, idx)

    if isempty(x) || isempty(idx) || idx > numel(x)
        y = NaN;
    else
        y = x(idx);
    end

end

function printTriggerConditions(C)

    fprintf('\nTest conditions at trigger:\n');
    fprintf('  Trigger time:        %.9f s\n', C.time_s);
    fprintf('  High-rate index:     %.0f\n', C.sampleIndex_highRate);
    fprintf('  Ps:                  %.6f kPa\n', C.Ps_kPa);
    fprintf('  Pt:                  %.6f kPa\n', C.Pt_kPa);
    fprintf('  T_inf:               %.6f K\n', C.Tinf_K);
    fprintf('  T0:                  %.6f K\n', C.T0_K);
    fprintf('  Mach:                %.6f\n', C.Mach);
    fprintf('  U_inf:               %.6f m/s\n', C.Uinf_mps);
    fprintf('  rho_inf:             %.9f kg/m^3\n', C.rhoinf_kgm3);

end

function plotIfAvailable(x, y, displayName, lineWidth, lineSpec)

    if nargin < 5
        lineSpec = '-';
    end

    if ~isempty(x) && ~isempty(y)
        plot(x, y, lineSpec, ...
            'LineWidth', lineWidth, ...
            'DisplayName', displayName);
    end

end

function addTriggerLine(triggerTime_s)

    if ~isempty(triggerTime_s)
        xline(triggerTime_s, '--', ...
            sprintf('Trigger = %.6f s', triggerTime_s), ...
            'LineWidth', 1.2, ...
            'LabelOrientation', 'horizontal', ...
            'LabelVerticalAlignment', 'bottom', ...
            'DisplayName', 'Trigger');
    end

end