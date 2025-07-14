% Enhanced DBS Spectrum Analysis Script with Lead Location and Hemisphere Detection
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

fprintf('Selected %d file(s) for analysis\n', length(fileNames));

%% 2) Load and concatenate all JSON files
all_BrainSenseTimeDomain = [];
all_BrainSenseLfp = [];
successful_files = {};
failed_files = {};

% Store lead location and hemisphere info extracted from files
file_lead_locations = {};
file_hemispheres = {};

for fileIdx = 1:length(fileNames)
    fprintf('Loading file %d/%d: %s\n', fileIdx, length(fileNames), fileNames{fileIdx});
    
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
        lead_location = 'Unknown';
        hemisphere = 'Unknown';
        
        % Stack-based recursive search through JSON structure
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
                                lead_location = 'GPi';
                            elseif contains(lead_str, 'stn', 'IgnoreCase', true)
                                lead_location = 'STN';
                            end
                        end
                    end
                    
                    % Check for Hemisphere
                    if strcmpi(fieldName, 'Hemisphere')
                        if ischar(fieldValue) || isstring(fieldValue)
                            hemisphere_str = lower(char(fieldValue));
                            if contains(hemisphere_str, 'left', 'IgnoreCase', true)
                                hemisphere = 'Left';
                            elseif contains(hemisphere_str, 'right', 'IgnoreCase', true)
                                hemisphere = 'Right';
                            end
                        end
                    end
                    
                    % Add nested structures to search stack
                    if isstruct(fieldValue)
                        search_stack{end+1} = fieldValue;
                    elseif iscell(fieldValue)
                        % Add cell array elements to search stack
                        for cellIdx = 1:numel(fieldValue)
                            if isstruct(fieldValue{cellIdx})
                                search_stack{end+1} = fieldValue{cellIdx};
                            end
                        end
                    end
                end
            elseif iscell(current_obj)
                % If current object is a cell array, add its struct elements
                for cellIdx = 1:numel(current_obj)
                    if isstruct(current_obj{cellIdx})
                        search_stack{end+1} = current_obj{cellIdx};
                    end
                end
            end
        end
        
        % Store the extracted information
        file_lead_locations{fileIdx} = lead_location;
        file_hemispheres{fileIdx} = hemisphere;
        
        fprintf('  Lead Location: %s, Hemisphere: %s\n', lead_location, hemisphere);
        
        % Extract the two required fields
        if ~isfield(data,'BrainSenseTimeDomain') || ~isfield(data,'BrainSenseLfp')
            warning('File %s missing required fields. Skipping.', fileNames{fileIdx});
            failed_files{end+1} = sprintf('%s (missing required fields)', fileNames{fileIdx});
            continue;
        end
        
        successful_files{end+1} = fileNames{fileIdx};
        
    catch ME
        % Detailed error reporting
        fprintf('\n*** ERROR loading %s ***\n', fileNames{fileIdx});
        fprintf('Error message: %s\n', ME.message);
        
        % Try to provide more specific error information
        if contains(ME.message, 'JSON syntax error')
            fprintf('This appears to be a malformed JSON file.\n');
            
            % Check file size
            fileInfo = dir(fullPath);
            fprintf('File size: %d bytes\n', fileInfo.bytes);
            
            % Show first and last few characters to help diagnose
            if length(jsonText) > 200
                fprintf('First 100 characters: %s\n', jsonText(1:min(100,end)));
                fprintf('Last 100 characters: %s\n', jsonText(max(1,end-99):end));
            else
                fprintf('File content: %s\n', jsonText);
            end
        end
        
        failed_files{end+1} = sprintf('%s (%s)', fileNames{fileIdx}, ME.message);
        fprintf('Skipping this file and continuing with others...\n\n');
        continue;
    end
    
    % Convert to table if needed
    if ~istable(data.BrainSenseTimeDomain)
        currentTD = struct2table(data.BrainSenseTimeDomain);
    else
        currentTD = data.BrainSenseTimeDomain;
    end
    
    % Add file source and metadata information
    currentTD.SourceFile = repmat({fileNames{fileIdx}}, height(currentTD), 1);
    currentTD.LeadLocation = repmat({lead_location}, height(currentTD), 1);
    currentTD.Hemisphere = repmat({hemisphere}, height(currentTD), 1);
    
    % Concatenate TimeDomain data
    if isempty(all_BrainSenseTimeDomain)
        all_BrainSenseTimeDomain = currentTD;
    else
        % Ensure matching columns
        commonVars = intersect(all_BrainSenseTimeDomain.Properties.VariableNames, ...
                               currentTD.Properties.VariableNames);
        all_BrainSenseTimeDomain = [all_BrainSenseTimeDomain(:, commonVars); ...
                                    currentTD(:, commonVars)];
    end
    
    % Concatenate LFP data
    if isempty(all_BrainSenseLfp)
        all_BrainSenseLfp = data.BrainSenseLfp(:);
    else
        all_BrainSenseLfp = [all_BrainSenseLfp; data.BrainSenseLfp(:)];
    end
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
tdTimes = datetime(BrainSenseTimeDomain.FirstPacketDateTime, ...
    'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
[~, sortIdx] = sort(tdTimes);
BrainSenseTimeDomain = BrainSenseTimeDomain(sortIdx, :);

% Sort LFP data as well
lfpTimes = datetime({BrainSenseLfp.FirstPacketDateTime}, ...
    'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
[~, sortIdx] = sort(lfpTimes);
BrainSenseLfp = BrainSenseLfp(sortIdx);

%% 3) Initialize StimRateHz
nTD = height(BrainSenseTimeDomain);
BrainSenseTimeDomain.StimRateHz = nan(nTD,1);

% Parse timestamps once
tdTimes = datetime(BrainSenseTimeDomain.FirstPacketDateTime, ...
    'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
nLFP = numel(BrainSenseLfp);
lfpTimes = datetime({BrainSenseLfp.FirstPacketDateTime}, ...
    'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC')';

% Loop & match ±1.5 s to get StimRateHz
tol = seconds(1.5);
for i = 1:nTD
    tdChan = BrainSenseTimeDomain.Channel{i};
    tdTime = tdTimes(i);
    
    % Find the first LFP row within tolerance
    dt = abs(lfpTimes - tdTime);
    idx = find(dt <= tol, 1);
    if isempty(idx)
        continue;
    end
    
    % Split LFP channels, pick the one matching tdChan
    lfpChs = strsplit(BrainSenseLfp(idx).Channel, ',');
    if ~any(strcmp(lfpChs, tdChan))
        continue;
    end
    
    % Grab the correct side's RateInHertz
    snap = BrainSenseLfp(idx).TherapySnapshot;
    if endsWith(tdChan,'_LEFT')
        rate = snap.Left.RateInHertz;
    else
        rate = snap.Right.RateInHertz;
    end
    
    BrainSenseTimeDomain.StimRateHz(i) = rate;
end

%% 4) Group data by Channel and StimRateHz (like original script)
fprintf('\nGrouping data by channel and stimulation rate...\n');

% Create unique combinations of Channel and StimRateHz
channels = BrainSenseTimeDomain.Channel;
stimRates = BrainSenseTimeDomain.StimRateHz;

% Create a combined identifier
channelRateCombos = {};
comboIndices = {};  % Store indices for each combo
comboLeadLocations = {};  % Store lead locations for each combo
comboHemispheres = {};  % Store hemispheres for each combo

for i = 1:nTD
    if ~isnan(stimRates(i))
        combo = sprintf('%s @ %.1f Hz', channels{i}, stimRates(i));
        
        % Check if this combo already exists
        existingIdx = find(strcmp(channelRateCombos, combo));
        if isempty(existingIdx)
            % New combo
            channelRateCombos{end+1} = combo;
            comboIndices{end+1} = i;
            
            % Get lead location and hemisphere for this combo
            lead_loc = BrainSenseTimeDomain.LeadLocation{i};
            hemisphere = BrainSenseTimeDomain.Hemisphere{i};
            comboLeadLocations{end+1} = lead_loc;
            comboHemispheres{end+1} = hemisphere;
        else
            % Add to existing combo
            comboIndices{existingIdx} = [comboIndices{existingIdx}, i];
        end
    end
end

% Display available combinations with file information
fprintf('\nAvailable Channel/StimRate combinations:\n');
for i = 1:length(channelRateCombos)
    % Get unique source files for this combination
    comboFiles = unique(BrainSenseTimeDomain.SourceFile(comboIndices{i}));
    
    % Display with lead location and hemisphere info
    fprintf('%d. %s | %s %s (%d segments from %d file(s))\n', i, ...
        channelRateCombos{i}, ...
        comboLeadLocations{i}, ...
        comboHemispheres{i}, ...
        length(comboIndices{i}), ...
        length(comboFiles));
end

% Prompt user for selection
selection = input('\nEnter the number of the combination you want to analyze: ');

if selection < 1 || selection > length(channelRateCombos)
    error('Invalid selection.');
end

selectedCombo = channelRateCombos{selection};
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
selected_stim_rate = stimRates(selectedIndices(1)); % Get stim rate from first selected index

% Calculate gamma band as 1/2 the stim rate (±1 Hz)
if ~isnan(selected_stim_rate)
    gamma_center = selected_stim_rate / 2;
    gamma_band = [gamma_center - 1, gamma_center + 1];
    fprintf('\nSelected stimulation rate: %.1f Hz\n', selected_stim_rate);
    fprintf('Calculated gamma band: [%.1f, %.1f] Hz (targeting %.1f Hz subharmonic)\n', ...
        gamma_band(1), gamma_band(2), gamma_center);
else
    error('No valid stimulation rate found for selected channel');
end

%% 6) Process epochs and compute spectra
% Parameters
epoch_duration_sec = 30;  % 30-second epochs
fs = 250;  % Sampling frequency (Hz)
welch_window_sec = 1;
welch_overlap = 0.5;
welch_window_samples = welch_window_sec * fs;
welch_overlap_samples = floor(welch_window_samples * welch_overlap);

fprintf('\n=== Processing Epochs for Spectrum Analysis ===\n');
fprintf('Epoch duration: %d seconds\n', epoch_duration_sec);

% Concatenate all time domain data for selected segments
total_samples = sum( cellfun(@(sig) numel(sig), ...
                     BrainSenseTimeDomain.TimeDomainData(selectedIndices)) );

% Preallocate numeric and datetime arrays
all_signals      = zeros(total_samples,1);
all_times        = NaT(total_samples,1,'TimeZone','UTC');   % <-- specify UTC here
all_file_sources = cell(total_samples,1);

% A running pointer into where we’ll write next
write_ptr = 1;

for idx = selectedIndices
    % Get the time domain data
    if ismember('TimeDomainData', BrainSenseTimeDomain.Properties.VariableNames) && ...
       ~isempty(BrainSenseTimeDomain.TimeDomainData{idx})
        
        signal = BrainSenseTimeDomain.TimeDomainData{idx};
        
        % Ensure signal is a column vector
        if size(signal, 2) > size(signal, 1)
            signal = signal';
        end
        
        % Get the start time for this segment
        startTime = tdTimes(idx);
        
        % Create time vector - each sample is 4ms apart
        nSamples = length(signal);
        sample_interval_ms = 4;
        
        % Generate timestamps for each sample
        sampleTimes = startTime + milliseconds((0:nSamples-1) * sample_interval_ms);
        
        % Track source file
        sourceFile = BrainSenseTimeDomain.SourceFile{idx};
        
        % Append to arrays
        nS = numel(signal);
        range = write_ptr:(write_ptr + nS - 1);

        all_signals    (range) = signal(:);
        all_times      (range) = sampleTimes(:);
        all_file_sources(range) = repmat({sourceFile},nS,1);

        write_ptr = write_ptr + nS;
    end
end

% Sort by time
[all_times, sort_idx] = sort(all_times);
all_signals = all_signals(sort_idx);
all_file_sources = all_file_sources(sort_idx);

% Check if we have data
if isempty(all_times) || isempty(all_signals)
    error('No time domain data available for the selected configuration.');
end

% Convert to relative time
start_time = all_times(1);
relative_time_sec = seconds(all_times - start_time);
total_duration_sec = relative_time_sec(end);
num_epochs = floor(total_duration_sec / epoch_duration_sec);

fprintf('Total recording duration: %.1f seconds\n', total_duration_sec);
fprintf('Number of potential epochs: %d\n', num_epochs);

%% 7) Compute spectra for all valid epochs
%–– build list of valid epochs ––
valid_epochs = false(num_epochs,1);
for e = 1:num_epochs
    % compute epoch start/end times exactly as before…
    epoch_start_time = start_time + seconds((e-1)*epoch_duration_sec);
    epoch_end_time   = epoch_start_time + seconds(epoch_duration_sec);
    mask = all_times>=epoch_start_time & all_times<epoch_end_time;
    sig = all_signals(mask);
    if isempty(sig), continue; end

    dt_ms = milliseconds(diff(all_times(mask)));
    dropped = any(dt_ms > (4 + 1));  % tolerance
    count_ok = abs(numel(sig) - epoch_duration_sec*fs) <= 0.01*epoch_duration_sec*fs;
    if ~dropped && count_ok
        valid_epochs(e) = true;
    end
end

% pre-allocate
idxs = find(valid_epochs);
nValid = numel(idxs);
[~, f] = pwelch(all_signals(1:fs*epoch_duration_sec), ...
                hamming(welch_window_samples), ...
                welch_overlap_samples, [], fs);
nFreq = numel(f);
all_psd_db = nan(nValid, nFreq);
valid_epoch_count = 0;

% loop only over valid ones
for vi = 1:nValid
    e = idxs(vi);
    epoch_start_time = start_time + seconds((e-1)*epoch_duration_sec);
    mask = all_times>=epoch_start_time & all_times<epoch_start_time+seconds(epoch_duration_sec);
    sig = all_signals(mask);
    [pxx,~] = pwelch(sig, hamming(welch_window_samples), ...
                    welch_overlap_samples, [], fs);
    all_psd_db(vi,:) = 10*log10(pxx)';
    valid_epoch_count = valid_epoch_count + 1;
    if mod(vi,10)==0
        fprintf('Processed valid epoch %d/%d\n', vi, nValid);
        drawnow;
    end
end

all_freqs = f;
fprintf('\nValid epochs processed: %d/%d\n', valid_epoch_count, num_epochs);


fprintf('\nValid epochs processed: %d/%d\n', valid_epoch_count, num_epochs);

%% 8) Plot spectrum with lead location and hemisphere in title
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
        valid_epoch_count, length(selectedFiles));
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
        strrep(selectedCombo, '_', '\_'), total_duration_sec/60, valid_epoch_count);
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
    spectrum_results.num_epochs = valid_epoch_count;
    spectrum_results.source_files = selectedFiles;
    
    assignin('base', 'spectrum_results', spectrum_results);
    
    fprintf('\n=== Analysis Complete ===\n');
    fprintf('Results saved to workspace as ''spectrum_results''\n');
    fprintf('Channel: %s\n', selectedCombo);
    fprintf('Lead Location: %s\n', selectedLeadLocation);
    fprintf('Hemisphere: %s\n', selectedHemisphere);
    fprintf('Stimulation Rate: %.1f Hz\n', selected_stim_rate);
    fprintf('Valid epochs: %d\n', valid_epoch_count);
    fprintf('Files processed: %d\n', length(selectedFiles));
    
else
    error('No valid epochs found for spectral analysis.');
end