% Optimized DBS Spectrum Analysis Script with Lead Location and Hemisphere Detection
% This script processes BrainSense JSON files, concatenates data by channel,
% and creates spectrum plots with lead location (STN/GPi) and hemisphere (Left/Right) info

clear; clc; close all;

%% 1) Prompt for multiple JSON files
[fileNames, filePath] = uigetfile('*.json','Select BrainSense JSON file(s)', 'MultiSelect', 'on');

% Handle single vs multiple file selection
if isequal(fileNames,0)
    error('No file selected. Exiting.');
end
% Convert single file to cell array for consistent handling
if ~iscell(fileNames)
    fileNames = {fileNames};
end
numFiles = length(fileNames);
fprintf('Selected %d file(s) for analysis\n', numFiles);

%% 2) Load and concatenate all JSON files - OPTIMIZED and CORRECTED
% Pre-allocate cell arrays to store data before concatenating into tables/structs
all_BrainSenseTimeDomain_cells = cell(numFiles, 1);
all_BrainSenseLfp_cells = cell(numFiles, 1);

successful_files = {};
failed_files = {};

for fileIdx = 1:numFiles
    fprintf('Loading file %d/%d: %s\n', fileIdx, numFiles, fileNames{fileIdx});

    current_lead_location = 'Unknown';
    current_hemisphere = 'Unknown';
    currentTD = []; % Initialize for current file's TimeDomain data
    currentLFP = []; % Initialize for current file's LFP data

    try
        % Read JSON with error handling
        fullPath = fullfile(filePath, fileNames{fileIdx});
        jsonText = fileread(fullPath);

        % Check if file is empty
        if isempty(strtrim(jsonText))
            error('File is empty');
        end

        % Try to decode JSON
        data = jsondecode(jsonText);

        % Extract lead location and hemisphere with recursive search
        search_stack = {data};

        while ~isempty(search_stack)
            current_obj = search_stack{end};
            search_stack(end) = [];

            if isstruct(current_obj)
                fields = fieldnames(current_obj);
                for idx = 1:length(fields)
                    fieldName = fields{idx};
                    fieldValue = current_obj.(fieldName);

                    % Check for LeadLocation
                    if strcmpi(fieldName, 'LeadLocation') || strcmpi(fieldName, 'Lead_Location')
                        if ischar(fieldValue) || isstring(fieldValue)
                            lead_str = lower(char(fieldValue));
                            if contains(lead_str, 'gpi', 'IgnoreCase', true)
                                current_lead_location = 'GPi';
                            elseif contains(lead_str, 'stn', 'IgnoreCase', true)
                                current_lead_location = 'STN';
                            end
                        end
                    end

                    % Check for Hemisphere
                    if strcmpi(fieldName, 'Hemisphere')
                        if ischar(fieldValue) || isstring(fieldValue)
                            hemisphere_str = lower(char(fieldValue));
                            if contains(hemisphere_str, 'left', 'IgnoreCase', true)
                                current_hemisphere = 'Left';
                            elseif contains(hemisphere_str, 'right', 'IgnoreCase', true)
                                current_hemisphere = 'Right';
                            end
                        end
                    end

                    % Add nested structures to search stack
                    if isstruct(fieldValue)
                        search_stack{end+1} = fieldValue;
                    elseif iscell(fieldValue)
                        for cellIdx = 1:numel(fieldValue)
                            if isstruct(fieldValue{cellIdx})
                                search_stack{end+1} = fieldValue{cellIdx};
                            end
                        end
                    end
                end
            elseif iscell(current_obj)
                for cellIdx = 1:numel(current_obj)
                    if isstruct(current_obj{cellIdx})
                        search_stack{end+1} = current_obj{cellIdx};
                    end
                end
            end
        end

        fprintf('  Lead Location: %s, Hemisphere: %s\n', current_lead_location, current_hemisphere);

        % Extract the two required fields
        if ~isfield(data,'BrainSenseTimeDomain') || ~isfield(data,'BrainSenseLfp')
            warning('File %s missing required fields. Skipping.', fileNames{fileIdx});
            failed_files{end+1} = sprintf('%s (missing required fields)', fileNames{fileIdx});
            continue;
        end

        successful_files{end+1} = fileNames{fileIdx};

        % Convert to table if needed
        if ~istable(data.BrainSenseTimeDomain)
            currentTD = struct2table(data.BrainSenseTimeDomain);
        else
            currentTD = data.BrainSenseTimeDomain;
        end

        currentLFP = data.BrainSenseLfp(:);

        % --- IMPORTANT FIX: Add metadata columns to currentTD here ---
        % These columns are essential for subsequent steps and should always be present.
        nRows = height(currentTD);
        currentTD.SourceFile = repmat({fileNames{fileIdx}}, nRows, 1); % Use cell array for strings
        currentTD.LeadLocation = repmat({current_lead_location}, nRows, 1);
        currentTD.Hemisphere = repmat({current_hemisphere}, nRows, 1);
        % ---------------------------------------------------------------

        % Store the data for this file in the pre-allocated cell arrays
        all_BrainSenseTimeDomain_cells{fileIdx} = currentTD;
        all_BrainSenseLfp_cells{fileIdx} = currentLFP;

    catch ME
        % Detailed error reporting
        fprintf('\n*** ERROR loading %s ***\n', fileNames{fileIdx});
        fprintf('Error message: %s\n', ME.message);

        if contains(ME.message, 'JSON syntax error')
            fprintf('This appears to be a malformed JSON file.\n');
            fileInfo = dir(fullPath);
            fprintf('File size: %d bytes\n', fileInfo.bytes);
            if length(jsonText) > 200
                fprintf('First 100 characters: %s\n', jsonText(1:min(100,end)));
                fprintf('Last 100 characters: %s\n', jsonText(max(1,end-99):end));
            else
                fprintf('File content: %s\n', jsonText);
            end
        end

        failed_files{end+1} = sprintf('%s (%s)', fileNames{fileIdx}, ME.message);
        fprintf('Skipping this file and continuing with others...\n\n');
        % Set cells to empty for failed files to avoid issues during final concatenation
        all_BrainSenseTimeDomain_cells{fileIdx} = [];
        all_BrainSenseLfp_cells{fileIdx} = [];
        continue; % Skip to next file
    end
end

% Filter out empty entries (from failed files) before final concatenation
all_BrainSenseTimeDomain_cells = all_BrainSenseTimeDomain_cells(~cellfun(@isempty, all_BrainSenseTimeDomain_cells));
all_BrainSenseLfp_cells = all_BrainSenseLfp_cells(~cellfun(@isempty, all_BrainSenseLfp_cells));


% Now, concatenate all collected data once
if ~isempty(all_BrainSenseTimeDomain_cells)
    % When concatenating tables with varying columns, it's safer to ensure
    % all tables have all *potential* columns, filling missing ones with
    % default values, or only concatenate common ones and then add specific
    % new columns.
    % The previous approach of finding commonVars might have excluded our new columns.
    % Let's use `vertcat` directly on the cell array, as the new columns
    % (`SourceFile`, `LeadLocation`, `Hemisphere`) are now guaranteed to be
    % added to each `currentTD` table before it's stored.
    all_BrainSenseTimeDomain = vertcat(all_BrainSenseTimeDomain_cells{:});
    all_BrainSenseLfp = vertcat(all_BrainSenseLfp_cells{:});
else
    all_BrainSenseTimeDomain = table(); % Empty table if no valid data
    all_BrainSenseLfp = struct([]); % Empty struct array
end


% Report loading summary
fprintf('\n=== File Loading Summary ===\n');
fprintf('Successfully loaded: %d file(s)\n', length(successful_files));
if ~isempty(successful_files)
    for i = 1:length(successful_files)
        fprintf('  ✓ %s\n', successful_files{i});
    end
end
if ~isempty(failed_files)
    fprintf('\nFailed to load: %d file(s)\n', length(failed_files));
    for i = 1:length(failed_files)
        fprintf('  ✗ %s\n', failed_files{i});
    end
end

% Check if we have any data to work with
if isempty(all_BrainSenseTimeDomain)
    error('No valid data was loaded from any file. Please check your JSON files.');
end

% Use concatenated data for the rest of the analysis
BrainSenseTimeDomain = all_BrainSenseTimeDomain;
BrainSenseLfp = all_BrainSenseLfp;

fprintf('\nTotal data loaded:\n');
fprintf('- TimeDomain entries: %d\n', height(BrainSenseTimeDomain));
fprintf('- LFP entries: %d\n', length(BrainSenseLfp));

% Sort data by timestamp to ensure chronological order
% Parse all timestamps at once for TimeDomain and LFP
tdTimes = datetime(BrainSenseTimeDomain.FirstPacketDateTime, ...
    'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
[~, sortIdxTD] = sort(tdTimes);
BrainSenseTimeDomain = BrainSenseTimeDomain(sortIdxTD, :);

% Check if BrainSenseLfp is not empty before attempting to sort
if ~isempty(BrainSenseLfp)
    lfpTimes = datetime({BrainSenseLfp.FirstPacketDateTime}, ...
        'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC')';
    [~, sortIdxLFP] = sort(lfpTimes);
    BrainSenseLfp = BrainSenseLfp(sortIdxLFP);
else
    lfpTimes = NaT(0,1,'TimeZone','UTC'); % Initialize empty lfpTimes if no LFP data
end

%% 3) Initialize StimRateHz and match - OPTIMIZED
nTD = height(BrainSenseTimeDomain);
BrainSenseTimeDomain.StimRateHz = nan(nTD,1);

if isempty(BrainSenseLfp)
    warning('No LFP data found. Cannot determine StimRateHz.');
else
    tol = seconds(1.5); % Tolerance for matching

    % Convert LFP channel strings to cell array for easier comparison
    lfpChannels = {BrainSenseLfp.Channel};
    lfpTherapySnapshots = {BrainSenseLfp.TherapySnapshot}; % Extract therapy snapshots

    % Use a nested loop or find for matching
    for i = 1:nTD
        tdChan = BrainSenseTimeDomain.Channel{i};
        tdTime = tdTimes(i); % Use the already parsed datetime

        % Find LFP indices within tolerance
        dt_lfp = abs(lfpTimes - tdTime);
        possibleLFPIndices = find(dt_lfp <= tol);

        if isempty(possibleLFPIndices)
            continue;
        end

        % Iterate through possible LFP matches
        for j = 1:length(possibleLFPIndices)
            lfpIdx = possibleLFPIndices(j);

            % Check if the LFP channel matches the TD channel
            lfpChs = strsplit(lfpChannels{lfpIdx}, ','); % strsplit is relatively fast
            if any(strcmp(lfpChs, tdChan))
                snap = lfpTherapySnapshots{lfpIdx}; % Access pre-extracted snapshot
                if isfield(snap, 'Left') && endsWith(tdChan,'_LEFT')
                    rate = snap.Left.RateInHertz;
                elseif isfield(snap, 'Right') && endsWith(tdChan,'_RIGHT')
                    rate = snap.Right.RateInHertz;
                else
                    % Handle cases where channel might not end with _LEFT/_RIGHT
                    % or therapy snapshot doesn't have Left/Right fields
                    continue; % Skip if cannot determine hemisphere from channel name
                end
                BrainSenseTimeDomain.StimRateHz(i) = rate;
                break; % Found a match, move to next TimeDomain entry
            end
        end
    end
end


%% 4) Group data by Channel and StimRateHz (like original script) - Minor Optimization
fprintf('\nGrouping data by channel and stimulation rate...\n');

% Create unique combinations of Channel and StimRateHz
channels = BrainSenseTimeDomain.Channel;
stimRates = BrainSenseTimeDomain.StimRateHz;

% Filter out NaN stimRates for valid combos upfront
validDataMask = ~isnan(stimRates);
validChannels = channels(validDataMask);
validStimRates = stimRates(validDataMask);
validIndices = find(validDataMask); % Original indices of valid data

% Create unique combo strings
comboStrings = arrayfun(@(ch, sr) sprintf('%s @ %.1f Hz', ch{1}, sr), ...
                        validChannels, validStimRates, 'UniformOutput', false);

[uniqueCombos, ~, ic] = unique(comboStrings);

comboIndices = cell(length(uniqueCombos), 1);
comboLeadLocations = cell(length(uniqueCombos), 1);
comboHemispheres = cell(length(uniqueCombos), 1);

for i = 1:length(uniqueCombos)
    currentComboMask = (ic == i);
    originalIdxsForCombo = validIndices(currentComboMask); % Get original indices
    comboIndices{i} = originalIdxsForCombo;

    % Lead Location and Hemisphere for this combo (take the first one, assuming consistent within combo)
    firstIdxInCombo = originalIdxsForCombo(1);
    % These fields are now guaranteed to be in BrainSenseTimeDomain table
    comboLeadLocations{i} = BrainSenseTimeDomain.LeadLocation{firstIdxInCombo};
    comboHemispheres{i} = BrainSenseTimeDomain.Hemisphere{firstIdxInCombo};
end

% Display available combinations with file information
fprintf('\nAvailable Channel/StimRate combinations:\n');
for i = 1:length(uniqueCombos)
    % Get unique source files for this combination
    comboFiles = unique(BrainSenseTimeDomain.SourceFile(comboIndices{i}));

    % Display with lead location and hemisphere info
    fprintf('%d. %s | %s %s (%d segments from %d file(s))\n', i, ...
        uniqueCombos{i}, ...
        comboLeadLocations{i}, ...
        comboHemispheres{i}, ...
        length(comboIndices{i}), ...
        length(comboFiles));
end

% Prompt user for selection
selection = input('\nEnter the number of the combination you want to analyze: ');
if selection < 1 || selection > length(uniqueCombos)
    error('Invalid selection.');
end

selectedCombo = uniqueCombos{selection};
selectedIndices = comboIndices{selection};
selectedLeadLocation = comboLeadLocations{selection};
selectedHemisphere = comboHemispheres{selection};

fprintf('\nSelected: %s | %s %s\n', selectedCombo, selectedLeadLocation, selectedHemisphere);
fprintf('Processing %d segments...\n', length(selectedIndices));

% Display which files contribute to this selection
selectedFiles = unique(BrainSenseTimeDomain.SourceFile(selectedIndices));
fprintf('Data from files: ');
for i = 1:length(selectedFiles)
    fprintf('%s', selectedFiles{i});
    if i < length(selectedFiles)
        fprintf(', ');
    end
end
fprintf('\n');

%% 5) Calculate gamma band
% Already efficient
selected_stim_rate = stimRates(selectedIndices(1)); % Get stim rate from first selected index
if isnan(selected_stim_rate)
    error('No valid stimulation rate found for selected channel combination.');
end

gamma_center = selected_stim_rate / 2;
gamma_band = [gamma_center - 1, gamma_center + 1];
fprintf('\nSelected stimulation rate: %.1f Hz\n', selected_stim_rate);
fprintf('Calculated gamma band: [%.1f, %.1f] Hz (targeting %.1f Hz subharmonic)\n', ...
    gamma_band(1), gamma_band(2), gamma_center);

%% 6) Process epochs and compute spectra - OPTIMIZED
% Parameters
epoch_duration_sec = 30;  % 30-second epochs
fs = 250;  % Sampling frequency (Hz)
welch_window_sec = 1;
welch_overlap = 0.5;
welch_window_samples = welch_window_sec * fs;
welch_overlap_samples = floor(welch_window_samples * welch_overlap);

fprintf('\n=== Processing Epochs for Spectrum Analysis ===\n');
fprintf('Epoch duration: %d seconds\n', epoch_duration_sec);

% Concatenate all time domain data for selected segments - OPTIMIZED
% Pre-calculate total samples needed to preallocate
total_samples_for_selected = sum( cellfun(@(sig) numel(sig), ...
                                   BrainSenseTimeDomain.TimeDomainData(selectedIndices)) );

all_signals_selected      = zeros(total_samples_for_selected,1);
all_times_selected        = NaT(total_samples_for_selected,1,'TimeZone','UTC');
all_file_sources_selected = cell(total_samples_for_selected,1);

write_ptr = 1;
for idx = selectedIndices' % Iterate over selected indices
    signal = BrainSenseTimeDomain.TimeDomainData{idx};
    if isempty(signal), continue; end % Skip empty signals

    % Ensure signal is a column vector
    if size(signal, 2) > size(signal, 1)
        signal = signal';
    end

    startTime = tdTimes(idx); % Use the already parsed datetime
    nSamples = length(signal);
    sample_interval_ms = 4;

    % Generate timestamps for each sample using vectorized operation
    sampleTimes = startTime + milliseconds((0:nSamples-1) * sample_interval_ms);

    sourceFile = BrainSenseTimeDomain.SourceFile{idx};

    % Append to arrays
    nS = numel(signal);
    range = write_ptr:(write_ptr + nS - 1);
    all_signals_selected(range) = signal(:);
    all_times_selected(range) = sampleTimes(:);
    all_file_sources_selected(range) = repmat({sourceFile},nS,1);
    write_ptr = write_ptr + nS;
end

% Trim pre-allocated arrays if they weren't completely filled
all_signals_selected = all_signals_selected(1:write_ptr-1);
all_times_selected = all_times_selected(1:write_ptr-1);
all_file_sources_selected = all_file_sources_selected(1:write_ptr-1);

% Sort by time (if not already perfectly sorted by previous concatenation)
[all_times_selected, sort_idx_selected] = sort(all_times_selected);
all_signals_selected = all_signals_selected(sort_idx_selected);
all_file_sources_selected = all_file_sources_selected(sort_idx_selected);

% Check if we have data
if isempty(all_times_selected) || isempty(all_signals_selected)
    error('No time domain data available for the selected configuration.');
end

% Convert to relative time
start_time = all_times_selected(1);
relative_time_sec = seconds(all_times_selected - start_time);
total_duration_sec = relative_time_sec(end);
num_epochs = floor(total_duration_sec / epoch_duration_sec);

fprintf('Total recording duration: %.1f seconds\n', total_duration_sec);
fprintf('Number of potential epochs: %d\n', num_epochs);

%% 7) Compute spectra for all valid epochs - OPTIMIZED
% Calculate expected number of samples per epoch
expected_samples_per_epoch = epoch_duration_sec * fs;
sample_interval_sec = 4e-3; % 4ms converted to seconds

% Find valid epochs in a more vectorized manner
% Create epoch start times
epoch_start_times_series = start_time + seconds((0:num_epochs-1)*epoch_duration_sec);

valid_epochs_mask = false(num_epochs, 1);
epoch_signal_cells = cell(num_epochs, 1);

fprintf('Identifying valid epochs...\n');
for e_idx = 1:num_epochs
    epoch_start_time = epoch_start_times_series(e_idx);
    epoch_end_time = epoch_start_time + seconds(epoch_duration_sec);

    % Find indices within this epoch
    mask = (all_times_selected >= epoch_start_time) & (all_times_selected < epoch_end_time);
    sig = all_signals_selected(mask);

    if isempty(sig), continue; end

    % Check for dropped samples/gaps by verifying continuous time
    epoch_time_diffs = diff(all_times_selected(mask));
    % Check if any time difference is significantly greater than the expected 4ms
    dropped = any(abs(seconds(epoch_time_diffs) - sample_interval_sec) > 1e-4); % 0.1ms tolerance

    % Check if the number of samples is approximately correct
    count_ok = abs(numel(sig) - expected_samples_per_epoch) <= 0.01 * expected_samples_per_epoch; % 1% tolerance

    if ~dropped && count_ok
        valid_epochs_mask(e_idx) = true;
        epoch_signal_cells{e_idx} = sig;
    end
end

% Remove empty cells
epoch_signal_cells = epoch_signal_cells(valid_epochs_mask);
nValid = numel(epoch_signal_cells);

if nValid == 0
    error('No valid epochs found for spectral analysis.');
end

% Pre-allocate PSD matrix
% Use the first valid epoch to determine frequency points
[~, f] = pwelch(epoch_signal_cells{1}, ...
                hamming(welch_window_samples), ...
                welch_overlap_samples, [], fs);
nFreq = numel(f);
all_psd_db = nan(nValid, nFreq);

fprintf('Computing spectra for %d valid epochs...\n', nValid);
for vi = 1:nValid
    sig = epoch_signal_cells{vi};
    [pxx,~] = pwelch(sig, hamming(welch_window_samples), ...
                    welch_overlap_samples, [], fs);
    all_psd_db(vi,:) = 10*log10(pxx)';
    if mod(vi,10)==0 || vi == nValid
        fprintf('Processed valid epoch %d/%d\n', vi, nValid);
        drawnow;
    end
end

all_freqs = f;
fprintf('\nValid epochs processed: %d/%d (total potential epochs)\n', nValid, num_epochs);


%% 8) Plot spectrum with lead location and hemisphere in title
% This section is largely fine, no significant performance bottlenecks here.
if ~isempty(all_psd_db)
    figure('Position', [100, 100, 1000, 700]);

    % Calculate mean and std of spectra
    mean_psd_db = mean(all_psd_db, 1);
    std_psd_db = std(all_psd_db, 0, 1);

    % Plot mean spectrum with error bounds
    hold on;

    % Plot individual spectra in light gray
    for i = 1:size(all_psd_db, 1)
        plot(all_freqs, all_psd_db(i, :), 'Color', [0.8, 0.8, 0.8], ...
             'LineWidth', 0.5, 'HandleVisibility', 'off');
    end

    % Plot mean spectrum
    plot(all_freqs, mean_psd_db, 'b-', 'LineWidth', 2.5);

    % Plot ±1 std as shaded area
    fill([all_freqs; flipud(all_freqs)], ...
         [mean_psd_db + std_psd_db, fliplr(mean_psd_db - std_psd_db)]', ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    % Highlight gamma band
    ylims = ylim;
    patch([gamma_band(1) gamma_band(2) gamma_band(2) gamma_band(1)], ...
          [ylims(1) ylims(1) ylims(2) ylims(2)], ...
          [1 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    % Add vertical lines for gamma band
    line([gamma_band(1) gamma_band(1)], ylims, 'Color', 'r', ...
         'LineStyle', '--', 'LineWidth', 1.5);
    line([gamma_band(2) gamma_band(2)], ylims, 'Color', 'r', ...
         'LineStyle', '--', 'LineWidth', 1.5);

    % Labels and title
    xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Power Spectral Density (dB/Hz)', 'FontSize', 14, 'FontWeight', 'bold');

    % Create informative title with lead location and hemisphere
    title_str = sprintf('%s %s Hemisphere - Power Spectrum\n%s | Stim Rate: %.1f Hz | %d epochs from %d file(s)', ...
        selectedLeadLocation, selectedHemisphere, ...
        strrep(selectedCombo, '_', '\_'), selected_stim_rate, ...
        nValid, length(selectedFiles)); % Use nValid for epoch count
    title(title_str, 'FontSize', 16, 'FontWeight', 'bold');

    % Set frequency limits
    xlim([0, min(100, max(all_freqs))]);

    % Grid
    grid on;
    grid minor;

    % Legend
    legend({'Mean Spectrum', '±1 STD', ...
            sprintf('Gamma Band (%.1f-%.1f Hz)', gamma_band(1), gamma_band(2))}, ...
           'Location', 'northeast', 'FontSize', 12);

    % Add text annotation with key information
    annotation_text = sprintf('Channel: %s\nTotal duration: %.1f min\nValid epochs: %d', ...
        strrep(selectedCombo, '_', '\_'), total_duration_sec/60, nValid); % Use nValid
    annotation('textbox', [0.02, 0.85, 0.25, 0.1], 'String', annotation_text, ...
               'BackgroundColor', 'white', 'EdgeColor', 'black', ...
               'FontSize', 10, 'FitBoxToText', 'on');

    hold off;

    % Save results to workspace
    spectrum_results = struct();
    spectrum_results.frequencies = all_freqs;
    spectrum_results.mean_psd_db = mean_psd_db;
    spectrum_results.std_psd_db = std_psd_db;
    spectrum_results.all_psd_db = all_psd_db;
    spectrum_results.gamma_band = gamma_band;
    spectrum_results.channel = selectedCombo;
    spectrum_results.lead_location = selectedLeadLocation;
    spectrum_results.hemisphere = selectedHemisphere;
    spectrum_results.stim_rate = selected_stim_rate;
    spectrum_results.num_epochs = nValid; % Use nValid
    spectrum_results.source_files = selectedFiles;

    assignin('base', 'spectrum_results', spectrum_results);

    fprintf('\n=== Analysis Complete ===\n');
    fprintf('Results saved to workspace as ''spectrum_results''\n');
    fprintf('Channel: %s\n', selectedCombo);
    fprintf('Lead Location: %s\n', selectedLeadLocation);
    fprintf('Hemisphere: %s\n', selectedHemisphere);
    fprintf('Stimulation Rate: %.1f Hz\n', selected_stim_rate);
    fprintf('Valid epochs: %d\n', nValid); % Use nValid
    fprintf('Files processed: %d\n', length(selectedFiles));

else
    error('No valid epochs found for spectral analysis.');
end