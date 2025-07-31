clear; clc; close all;

% 1) File selection
[fileNames, filePath] = uigetfile('*.json','Select BrainSense JSON file(s)', 'MultiSelect', 'on');

if isequal(fileNames,0)
    error('No file selected. Exiting.');
end

if ~iscell(fileNames)
    fileNames = {fileNames};
end

fprintf('Selected %d file(s) for analysis\n', length(fileNames));

% 2) REVOLUTIONARY: Extract JSON sections without parsing entire files
fprintf('\n=== Revolutionary Selective JSON Section Extraction ===\n');

% Pre-check file sizes
file_info = [];
for i = 1:length(fileNames)
    fullPath = fullfile(filePath, fileNames{i});
    info = dir(fullPath);
    if isempty(info)
        warning('File %s not found', fileNames{i});
        continue;
    end
    file_info(end+1).name = fileNames{i};
    file_info(end).path = fullPath;
    file_info(end).size_mb = info.bytes / (1024^2);
end

fprintf('File sizes: ');
for i = 1:length(file_info)
    fprintf('%s (%.1fMB) ', file_info(i).name, file_info(i).size_mb);
end
fprintf('\n');

% Initialize storage
td_cells = cell(length(file_info), 1);
lfp_cells = cell(length(file_info), 1);
metadata_array = repmat(struct('success', false, 'filename', '', 'lead_location', 'Unknown', 'hemisphere', 'Unknown'), length(file_info), 1);

% Process each file with revolutionary extraction
for fileIdx = 1:length(file_info)
    fprintf('\nProcessing file %d/%d: %s (%.1fMB)\n', fileIdx, length(file_info), file_info(fileIdx).name, file_info(fileIdx).size_mb);
    
    % Initialize outputs
    metadata_temp = struct('success', false, 'filename', file_info(fileIdx).name, ...
                         'lead_location', 'Unknown', 'hemisphere', 'Unknown');
    td_table_temp = table();
    lfp_array_temp = [];
    
    try
        % Skip massive files
        if file_info(fileIdx).size_mb > 500
            warning('Skipping large file %s (%.1fMB)', file_info(fileIdx).name, file_info(fileIdx).size_mb);
            td_cells{fileIdx} = td_table_temp;
            lfp_cells{fileIdx} = lfp_array_temp;
            metadata_array(fileIdx) = metadata_temp;
            continue;
        end
        
        % BREAKTHROUGH: Revolutionary selective extraction
        fprintf('  Phase 1: Reading file as text...\n');
        tic;
        jsonText = fileread(file_info(fileIdx).path);
        fprintf('    File read in %.2f seconds\n', toc);
        
        if isempty(strtrim(jsonText))
            warning('File %s is empty', file_info(fileIdx).name);
            td_cells{fileIdx} = td_table_temp;
            lfp_cells{fileIdx} = lfp_array_temp;
            metadata_array(fileIdx) = metadata_temp;
            continue;
        end
        
        fprintf('  Phase 2: Fast metadata extraction...\n');
        tic;
        [lead_location, hemisphere] = extract_metadata_from_text(jsonText);
        metadata_temp.lead_location = lead_location;
        metadata_temp.hemisphere = hemisphere;
        fprintf('    Metadata extracted in %.2f seconds\n', toc);
        
        fprintf('  Phase 3: Extracting TimeDomain section...\n');
        tic;
        td_table_temp = extract_timedomain_section(jsonText, file_info(fileIdx).name, lead_location, hemisphere);
        fprintf('    TimeDomain extracted in %.2f seconds (%d entries)\n', toc, height(td_table_temp));
        
        fprintf('  Phase 4: Extracting LFP section...\n');
        tic;
        lfp_array_temp = extract_lfp_section(jsonText);
        fprintf('    LFP extracted in %.2f seconds (%d entries)\n', toc, length(lfp_array_temp));
        
        % Clear the massive jsonText immediately
        clear jsonText;
        
        if ~isempty(td_table_temp) && ~isempty(lfp_array_temp)
            metadata_temp.success = true;
            fprintf('  ✓ SUCCESS: Lead: %s, Hemisphere: %s\n', lead_location, hemisphere);
        else
            fprintf('  ✗ FAILED: Missing required data sections\n');
        end
        
    catch ME
        fprintf('  ✗ ERROR: %s\n', ME.message);
    end
    
    td_cells{fileIdx} = td_table_temp;
    lfp_cells{fileIdx} = lfp_array_temp;
    metadata_array(fileIdx) = metadata_temp;
end

% Filter successful files
success_mask = [metadata_array.success];
td_cells = td_cells(success_mask);
lfp_cells = lfp_cells(success_mask);
metadata_array = metadata_array(success_mask);
successful_count = sum(success_mask);

fprintf('\n=== Loading Summary ===\n');
fprintf('Successfully processed: %d/%d files\n', successful_count, length(file_info));

if successful_count == 0
    error('No valid data was loaded from any file. Please check your JSON files.');
end

% 3) Data concatenation with table structure standardization
fprintf('\n=== Data Concatenation with Structure Standardization ===\n');
tic;

% Filter non-empty data
non_empty_td = ~cellfun(@(x) isempty(x) || height(x) == 0, td_cells);
non_empty_lfp = ~cellfun(@isempty, lfp_cells);

if any(non_empty_td)
    fprintf('Standardizing TimeDomain table structures...\n');
    
    % Get all non-empty tables
    valid_td_tables = td_cells(non_empty_td);
    
    % Debug: Show table structures
    fprintf('Analyzing table structures across %d files:\n', length(valid_td_tables));
    for i = 1:length(valid_td_tables)
        fprintf('  File %d: %d rows, %d columns\n', i, height(valid_td_tables{i}), width(valid_td_tables{i}));
    end
    
    % Find all unique column names across all tables
    all_column_names = {};
    for i = 1:length(valid_td_tables)
        table_columns = valid_td_tables{i}.Properties.VariableNames;
        all_column_names = [all_column_names, table_columns];
    end
    unique_columns = unique(all_column_names);
    
    % Find common columns across ALL tables
    common_columns = valid_td_tables{1}.Properties.VariableNames;
    for i = 2:length(valid_td_tables)
        current_columns = valid_td_tables{i}.Properties.VariableNames;
        common_columns = intersect(common_columns, current_columns);
    end
    
    fprintf('Found %d unique columns, %d common columns across all tables\n', ...
        length(unique_columns), length(common_columns));
    
    % Try standardization first, fall back to common columns if it fails
    try
        % Standardize all tables to have the same columns
        standardized_tables = cell(length(valid_td_tables), 1);
        
        for i = 1:length(valid_td_tables)
            current_table = valid_td_tables{i};
            current_columns = current_table.Properties.VariableNames;
            
            % Find missing columns
            missing_columns = setdiff(unique_columns, current_columns);
            
            if ~isempty(missing_columns)
                fprintf('  Table %d: Adding %d missing columns: %s\n', i, length(missing_columns), strjoin(missing_columns, ', '));
                
                % Add missing columns with appropriate default values
                for j = 1:length(missing_columns)
                    col_name = missing_columns{j};
                    
                    % Determine appropriate default value based on column name and existing data
                    if contains(col_name, 'DateTime') || contains(col_name, 'Time')
                        current_table.(col_name) = NaT(height(current_table), 1, 'TimeZone', 'UTC');
                    elseif contains(col_name, 'Hz') || contains(col_name, 'Rate') || contains(col_name, 'Amplitude')
                        current_table.(col_name) = nan(height(current_table), 1);
                    elseif any(contains(current_columns, {'SourceFile', 'Channel', 'LeadLocation', 'Hemisphere'}))
                        % String/cell columns
                        current_table.(col_name) = repmat({''}, height(current_table), 1);
                    else
                        % Try to infer from existing columns
                        sample_col = current_table{:,1};
                        if iscell(sample_col) || isstring(sample_col)
                            current_table.(col_name) = repmat({''}, height(current_table), 1);
                        else
                            current_table.(col_name) = nan(height(current_table), 1);
                        end
                    end
                end
            end
            
            % Reorder columns to match the standard order
            current_table = current_table(:, unique_columns);
            standardized_tables{i} = current_table;
        end
        
        % Now concatenate the standardized tables
        fprintf('Concatenating standardized tables...\n');
        all_BrainSenseTimeDomain = vertcat(standardized_tables{:});
        
    catch ME
        fprintf('Standardization failed: %s\n', ME.message);
        fprintf('Falling back to common columns only...\n');
        
        % Fallback: use only common columns
        if isempty(common_columns)
            error('No common columns found across TimeDomain tables');
        end
        
        fprintf('Using %d common columns: %s\n', length(common_columns), strjoin(common_columns, ', '));
        
        % Extract only common columns from each table
        common_tables = cell(length(valid_td_tables), 1);
        for i = 1:length(valid_td_tables)
            common_tables{i} = valid_td_tables{i}(:, common_columns);
        end
        
        all_BrainSenseTimeDomain = vertcat(common_tables{:});
        clear common_tables;
    end
    
    fprintf('TimeDomain concatenated: %d total entries with %d columns\n', ...
        height(all_BrainSenseTimeDomain), width(all_BrainSenseTimeDomain));
else
    error('No TimeDomain data found');
end

if any(non_empty_lfp)
    fprintf('Standardizing LFP data structures...\n');
    
    % Get all non-empty LFP arrays
    valid_lfp_arrays = lfp_cells(non_empty_lfp);
    
    % Check if LFP structures are consistent
    first_lfp_fields = [];
    if ~isempty(valid_lfp_arrays{1})
        first_lfp_fields = fieldnames(valid_lfp_arrays{1}(1));
    end
    
    % Check for field consistency across all LFP arrays
    fields_consistent = true;
    for i = 2:length(valid_lfp_arrays)
        if ~isempty(valid_lfp_arrays{i})
            current_fields = fieldnames(valid_lfp_arrays{i}(1));
            if ~isequal(first_lfp_fields, current_fields)
                fields_consistent = false;
                fprintf('  Warning: LFP field structure differs between files\n');
                break;
            end
        end
    end
    
    if fields_consistent
        % Simple concatenation if structures are consistent
        all_BrainSenseLfp = vertcat(valid_lfp_arrays{:});
    else
        % Handle inconsistent structures by finding common fields
        fprintf('  Handling inconsistent LFP structures...\n');
        
        % Find all unique field names
        all_lfp_fields = {};
        for i = 1:length(valid_lfp_arrays)
            if ~isempty(valid_lfp_arrays{i})
                current_fields = fieldnames(valid_lfp_arrays{i}(1));
                all_lfp_fields = [all_lfp_fields; current_fields];
            end
        end
        unique_lfp_fields = unique(all_lfp_fields);
        
        % Standardize each LFP array
        standardized_lfp_arrays = cell(length(valid_lfp_arrays), 1);
        
        for i = 1:length(valid_lfp_arrays)
            if ~isempty(valid_lfp_arrays{i})
                current_array = valid_lfp_arrays{i};
                current_fields = fieldnames(current_array(1));
                missing_fields = setdiff(unique_lfp_fields, current_fields);
                
                % Add missing fields with default values
                if ~isempty(missing_fields)
                    for j = 1:length(current_array)
                        for k = 1:length(missing_fields)
                            field_name = missing_fields{k};
                            if contains(field_name, 'DateTime') || contains(field_name, 'Time')
                                current_array(j).(field_name) = NaT('TimeZone', 'UTC');
                            elseif contains(field_name, 'Hz') || contains(field_name, 'Rate')
                                current_array(j).(field_name) = NaN;
                            else
                                current_array(j).(field_name) = [];
                            end
                        end
                    end
                end
                standardized_lfp_arrays{i} = current_array;
            end
        end
        
        % Now concatenate the standardized arrays
        all_BrainSenseLfp = vertcat(standardized_lfp_arrays{:});
    end
    
    fprintf('LFP concatenated: %d total entries\n', length(all_BrainSenseLfp));
else
    error('No LFP data found');
end

clear td_cells lfp_cells valid_td_tables standardized_tables valid_lfp_arrays standardized_lfp_arrays;
fprintf('Data concatenation completed in %.2f seconds\n', toc);

% 4) StimRate matching
fprintf('\n=== StimRate Matching ===\n');
tic;

td_times = all_BrainSenseTimeDomain.DateTime;
lfp_times = [all_BrainSenseLfp.DateTime]';

td_times_num = datenum(td_times);
lfp_times_num = datenum(lfp_times);

fprintf('Matching %d TD entries with %d LFP entries...\n', length(td_times_num), length(lfp_times_num));

distance_threshold = 1.5 / (24 * 3600);

% Optimized nearest neighbor search
if length(td_times_num) * length(lfp_times_num) < 100000
    distance_matrix = abs(td_times_num(:) - lfp_times_num(:)');
    [distances, nearest_indices] = min(distance_matrix, [], 2);
    fprintf('Used vectorized distance matrix\n');
else
    % Chunked processing for large datasets
    chunk_size = 1000;
    nearest_indices = zeros(length(td_times_num), 1);
    distances = zeros(length(td_times_num), 1);
    
    for start_idx = 1:chunk_size:length(td_times_num)
        end_idx = min(start_idx + chunk_size - 1, length(td_times_num));
        chunk_indices = start_idx:end_idx;
        
        chunk_td_times = td_times_num(chunk_indices);
        distance_matrix_chunk = abs(chunk_td_times(:) - lfp_times_num(:)');
        [chunk_distances, chunk_nearest] = min(distance_matrix_chunk, [], 2);
        
        distances(chunk_indices) = chunk_distances;
        nearest_indices(chunk_indices) = chunk_nearest;
    end
    fprintf('Used chunked processing\n');
end

% Extract stimulation rates
all_BrainSenseTimeDomain.StimRateHz = nan(height(all_BrainSenseTimeDomain), 1);

valid_matches = distances <= distance_threshold;
valid_td_indices = find(valid_matches);
valid_lfp_indices = nearest_indices(valid_matches);

for i = 1:length(valid_td_indices)
    td_idx = valid_td_indices(i);
    lfp_idx = valid_lfp_indices(i);
    
    tdChan = all_BrainSenseTimeDomain.Channel{td_idx};
    lfpChs = strsplit(all_BrainSenseLfp(lfp_idx).Channel, ',');
    
    if any(strcmp(lfpChs, tdChan))
        snap = all_BrainSenseLfp(lfp_idx).TherapySnapshot;
        if endsWith(tdChan, '_LEFT')
            rate = snap.Left.RateInHertz;
        else
            rate = snap.Right.RateInHertz;
        end
        all_BrainSenseTimeDomain.StimRateHz(td_idx) = rate;
    end
end

fprintf('StimRate matching completed in %.2f seconds\n', toc);
fprintf('Found rates for %d/%d entries\n', sum(~isnan(all_BrainSenseTimeDomain.StimRateHz)), height(all_BrainSenseTimeDomain));

% 5) Channel analysis - Process ALL combinations automatically
fprintf('\n=== Channel Analysis - Processing ALL Combinations ===\n');

valid_mask = ~isnan(all_BrainSenseTimeDomain.StimRateHz);
valid_data = all_BrainSenseTimeDomain(valid_mask, :);

if isempty(valid_data)
    error('No valid stimulation rate data found.');
end

combo_strings = arrayfun(@(i) sprintf('%s @ %.1f Hz', ...
    valid_data.Channel{i}, valid_data.StimRateHz(i)), ...
    1:height(valid_data), 'UniformOutput', false);

[unique_combos, ~, combo_indices] = unique(combo_strings);

channelRateCombos = unique_combos;
comboIndices = cell(length(unique_combos), 1);
comboLeadLocations = cell(length(unique_combos), 1);
comboHemispheres = cell(length(unique_combos), 1);
comboStimRates = zeros(length(unique_combos), 1);

for i = 1:length(unique_combos)
    mask = combo_indices == i;
    indices_in_valid = find(mask);
    
    original_indices = find(valid_mask);
    comboIndices{i} = original_indices(indices_in_valid);
    
    first_idx = comboIndices{i}(1);
    comboLeadLocations{i} = all_BrainSenseTimeDomain.LeadLocation{first_idx};
    comboHemispheres{i} = all_BrainSenseTimeDomain.Hemisphere{first_idx};
    comboStimRates(i) = all_BrainSenseTimeDomain.StimRateHz(first_idx);
end

fprintf('\nFound %d Channel/StimRate combinations:\n', length(channelRateCombos));
for i = 1:length(channelRateCombos)
    comboFiles = unique(all_BrainSenseTimeDomain.SourceFile(comboIndices{i}));
    
    fprintf('%d. %s | %s %s (%d segments from %d file(s))\n', i, ...
        channelRateCombos{i}, ...
        comboLeadLocations{i}, ...
        comboHemispheres{i}, ...
        length(comboIndices{i}), ...
        length(comboFiles));
end

%% 6) Create subplot layout for all combinations
num_combos = length(channelRateCombos);

% Calculate optimal subplot layout
subplot_cols = ceil(sqrt(num_combos));
subplot_rows = ceil(num_combos / subplot_cols);

% Create single large figure
figure('Position', [50, 50, 300*subplot_cols, 250*subplot_rows], 'Name', 'DBS Spectrum Analysis - All Channels');

% Store all results
all_spectrum_results = cell(num_combos, 1);

% 7) Process each combination automatically
fprintf('\n=== Processing All %d Combinations ===\n', num_combos);

for combo_idx = 1:num_combos
    fprintf('\nProcessing combination %d/%d: %s\n', combo_idx, num_combos, channelRateCombos{combo_idx});
    
    selectedCombo = channelRateCombos{combo_idx};
    selectedIndices = comboIndices{combo_idx};
    selectedLeadLocation = comboLeadLocations{combo_idx};
    selectedHemisphere = comboHemispheres{combo_idx};
    selected_stim_rate = comboStimRates(combo_idx);
    
    % Calculate gamma band
    gamma_center = selected_stim_rate / 2;
    gamma_band = [gamma_center - 1, gamma_center + 1];
    
    % Get selected data
    selected_data = all_BrainSenseTimeDomain(selectedIndices, :);
    signal_cells = selected_data.TimeDomainData;
    time_cells = selected_data.DateTime;
    
    non_empty_mask = ~cellfun(@isempty, signal_cells);
    signal_cells = signal_cells(non_empty_mask);
    time_cells = time_cells(non_empty_mask);
    
    if isempty(signal_cells)
        fprintf('  WARNING: No valid signal data for %s, skipping...\n', selectedCombo);
        continue;
    end
    
    % Calculate segment durations
    segment_durations = zeros(length(signal_cells), 1);
    fs = 250; % Sampling frequency
    sample_interval_ms = 4; % 4ms = 250Hz
    
    for i = 1:length(signal_cells)
        signal = signal_cells{i};
        if size(signal, 2) > size(signal, 1)
            signal = signal';
        end
        nSamples = length(signal);
        duration_sec = nSamples * sample_interval_ms / 1000;
        segment_durations(i) = duration_sec;
    end
    
    total_actual_duration = sum(segment_durations);
    
    % Process epochs
    epoch_duration_sec = 10;
    welch_window_sec = 1;
    welch_overlap = 0.5;
    welch_window_samples = welch_window_sec * fs;
    welch_overlap_samples = floor(welch_window_samples * welch_overlap);
    
    % Collect all valid epochs from all segments
    all_epoch_signals = {};
    total_epochs = 0;
    
    for seg_idx = 1:length(signal_cells)
        signal = signal_cells{seg_idx};
        if size(signal, 2) > size(signal, 1)
            signal = signal';
        end
        
        segment_duration_sec = length(signal) * sample_interval_ms / 1000;
        num_epochs_this_segment = floor(segment_duration_sec / epoch_duration_sec);
        
        % Extract epochs from this segment
        for e = 1:num_epochs_this_segment
            start_sample = (e-1) * epoch_duration_sec * fs + 1;
            end_sample = min(start_sample + epoch_duration_sec * fs - 1, length(signal));
            
            if (end_sample - start_sample + 1) >= 0.99 * epoch_duration_sec * fs
                epoch_signal = signal(start_sample:end_sample);
                all_epoch_signals{end+1} = epoch_signal;
                total_epochs = total_epochs + 1;
            end
        end
    end
    
    fprintf('  Found %d epochs (%.1f min total)\n', total_epochs, total_actual_duration/60);
    
    if total_epochs == 0
        fprintf('  WARNING: No valid epochs for %s, skipping...\n', selectedCombo);
        continue;
    end
    
    % Calculate spectra for all epochs
    [~, f] = pwelch(all_epoch_signals{1}, ...
                    hamming(welch_window_samples), ...
                    welch_overlap_samples, [], fs);
    
    all_psd_db = zeros(total_epochs, length(f));
    
    for i = 1:total_epochs
        sig = all_epoch_signals{i};
        [pxx, ~] = pwelch(sig, hamming(welch_window_samples), ...
                          welch_overlap_samples, [], fs);
        all_psd_db(i, :) = 10 * log10(pxx)';
    end
    
    % Create subplot for this combination
    subplot(subplot_rows, subplot_cols, combo_idx);
    
    mean_psd_db = mean(all_psd_db, 1);
    std_psd_db = std(all_psd_db, 0, 1);
    
    hold on;
    
    % Plot individual epoch spectra as light gray lines (reduced opacity for cleaner look)
    for i = 1:min(20, total_epochs) % Limit to 20 traces for cleaner plots
        plot(f, all_psd_db(i, :), 'Color', [0.9, 0.9, 0.9], ...
             'LineWidth', 0.2, 'HandleVisibility', 'off');
    end
    
    % Plot mean spectrum
    plot(f, mean_psd_db, 'b-', 'LineWidth', 1.5);
    
    % Plot confidence interval
    fill([f; flipud(f)], ...
         [mean_psd_db + std_psd_db, fliplr(mean_psd_db - std_psd_db)]', ...
         'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    
    % Highlight gamma band
    ylims = ylim;
    patch([gamma_band(1) gamma_band(2) gamma_band(2) gamma_band(1)], ...
          [ylims(1) ylims(1) ylims(2) ylims(2)], ...
          [1 0.8 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    line([gamma_band(1) gamma_band(1)], ylims, 'Color', 'r', ...
         'LineStyle', '--', 'LineWidth', 1);
    line([gamma_band(2) gamma_band(2)], ylims, 'Color', 'r', ...
         'LineStyle', '--', 'LineWidth', 1);
    
    xlabel('Frequency (Hz)', 'FontSize', 10);
    ylabel('PSD (dB/Hz)', 'FontSize', 10);
    
    title_str = sprintf('%s %s\n%s | %d epochs', ...
        selectedLeadLocation, selectedHemisphere, ...
        strrep(selectedCombo, '_', '\_'), total_epochs);
    title(title_str, 'FontSize', 10, 'FontWeight', 'bold');
    
    % MODIFIED: Set consistent x-axis limits and ticks for all subplots
    xlim([0, 100]);
    xticks(0:10:100);
    grid on;
    
    hold off;
    
    % Store results
    spectrum_result = struct();
    spectrum_result.frequencies = f;
    spectrum_result.mean_psd_db = mean_psd_db;
    spectrum_result.std_psd_db = std_psd_db;
    spectrum_result.all_psd_db = all_psd_db;
    spectrum_result.gamma_band = gamma_band;
    spectrum_result.channel = selectedCombo;
    spectrum_result.lead_location = selectedLeadLocation;
    spectrum_result.hemisphere = selectedHemisphere;
    spectrum_result.stim_rate = selected_stim_rate;
    spectrum_result.num_epochs = total_epochs;
    spectrum_result.total_duration_minutes = total_actual_duration/60;
    spectrum_result.epoch_duration_sec = epoch_duration_sec;
    
    all_spectrum_results{combo_idx} = spectrum_result;
end

% Add overall figure title
sgtitle('DBS Power Spectrum Analysis - All Channel/Rate Combinations', 'FontSize', 16, 'FontWeight', 'bold');

% Save all results
assignin('base', 'all_spectrum_results', all_spectrum_results);
assignin('base', 'channelRateCombos', channelRateCombos);

fprintf('\n=== COMPLETE ===\n');
fprintf('Processed %d channel combinations\n', num_combos);
fprintf('Results saved as ''all_spectrum_results'' and ''channelRateCombos''\n');

% EXTRACTION FUNCTIONS

function [lead_location, hemisphere] = extract_metadata_from_text(jsonText)
    % Ultra-fast metadata extraction using string search
    lead_location = 'Unknown';
    hemisphere = 'Unknown';
    
    jsonText_lower = lower(jsonText);
    
    % Lead location
    if contains(jsonText_lower, '"gpi"') || contains(jsonText_lower, 'leadlocationdef.gpi')
        lead_location = 'GPi';
    elseif contains(jsonText_lower, '"stn"') || contains(jsonText_lower, 'leadlocationdef.stn')
        lead_location = 'STN';
    end
    
    % Hemisphere
    if contains(jsonText_lower, '"left"') || contains(jsonText_lower, 'hemispherelocationdef.left')
        hemisphere = 'Left';
    elseif contains(jsonText_lower, '"right"') || contains(jsonText_lower, 'hemispherelocationdef.right')
        hemisphere = 'Right';
    end
end

function td_table = extract_timedomain_section(jsonText, filename, lead_location, hemisphere)
    % Extract only BrainSenseTimeDomain section and parse it
    td_table = table();
    
    try
        % Find BrainSenseTimeDomain section boundaries
        td_start_pattern = '"BrainSenseTimeDomain"';
        td_start = strfind(jsonText, td_start_pattern);
        
        if isempty(td_start)
            warning('BrainSenseTimeDomain not found in %s', filename);
            return;
        end
        
        % Find the start of the array after the field name
        colon_pos = strfind(jsonText(td_start(1):end), ':');
        if isempty(colon_pos)
            return;
        end
        
        array_start = td_start(1) + colon_pos(1);
        
        % Find matching bracket to get the complete array
        bracket_count = 0;
        in_string = false;
        escape_next = false;
        
        for i = array_start:length(jsonText)
            char = jsonText(i);
            
            if escape_next
                escape_next = false;
                continue;
            end
            
            if char == '\'
                escape_next = true;
                continue;
            end
            
            if char == '"'
                in_string = ~in_string;
                continue;
            end
            
            if ~in_string
                if char == '['
                    bracket_count = bracket_count + 1;
                elseif char == ']'
                    bracket_count = bracket_count - 1;
                    if bracket_count == 0
                        % Found the end of BrainSenseTimeDomain array
                        td_json_section = jsonText(array_start:i);
                        break;
                    end
                end
            end
        end
        
        if bracket_count ~= 0
            warning('Could not find complete BrainSenseTimeDomain section in %s', filename);
            return;
        end
        
        % Parse only this small section
        td_data = jsondecode(td_json_section);
        
        if ~isempty(td_data)
            if ~istable(td_data)
                td_table = struct2table(td_data);
            else
                td_table = td_data;
            end
            
            % Add metadata
            n_rows = height(td_table);
            td_table.SourceFile = repmat({filename}, n_rows, 1);
            td_table.LeadLocation = repmat({lead_location}, n_rows, 1);
            td_table.Hemisphere = repmat({hemisphere}, n_rows, 1);
            
            % Convert timestamps
            if iscell(td_table.FirstPacketDateTime)
                td_datetime_cells = td_table.FirstPacketDateTime;
            else
                td_datetime_cells = cellstr(string(td_table.FirstPacketDateTime));
            end
            
            td_table.DateTime = datetime(td_datetime_cells, ...
                'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
            
            % Sort by timestamp
            [~, sort_idx] = sort(td_table.DateTime);
            td_table = td_table(sort_idx, :);
        end
        
    catch ME
        warning('Error extracting TimeDomain from %s: %s', filename, ME.message);
    end
end

function lfp_array = extract_lfp_section(jsonText)
    % Extract only BrainSenseLfp section and parse it
    lfp_array = [];
    
    try
        % Find BrainSenseLfp section boundaries
        lfp_start_pattern = '"BrainSenseLfp"';
        lfp_start = strfind(jsonText, lfp_start_pattern);
        
        if isempty(lfp_start)
            warning('BrainSenseLfp not found');
            return;
        end
        
        % Find the start of the array
        colon_pos = strfind(jsonText(lfp_start(1):end), ':');
        if isempty(colon_pos)
            return;
        end
        
        array_start = lfp_start(1) + colon_pos(1);
        
        % Find matching bracket
        bracket_count = 0;
        in_string = false;
        escape_next = false;
        
        for i = array_start:length(jsonText)
            char = jsonText(i);
            
            if escape_next
                escape_next = false;
                continue;
            end
            
            if char == '\'
                escape_next = true;
                continue;
            end
            
            if char == '"'
                in_string = ~in_string;
                continue;
            end
            
            if ~in_string
                if char == '['
                    bracket_count = bracket_count + 1;
                elseif char == ']'
                    bracket_count = bracket_count - 1;
                    if bracket_count == 0
                        lfp_json_section = jsonText(array_start:i);
                        break;
                    end
                end
            end
        end
        
        if bracket_count ~= 0
            warning('Could not find complete BrainSenseLfp section');
            return;
        end
        
        % Parse only this small section
        lfp_data = jsondecode(lfp_json_section);
        
        if ~isempty(lfp_data)
            lfp_array = lfp_data(:);
            
            % Convert timestamps
            for i = 1:length(lfp_array)
                lfp_array(i).DateTime = datetime(lfp_array(i).FirstPacketDateTime, ...
                    'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
            end
            
            % Sort by timestamp
            lfp_times = [lfp_array.DateTime];
            [~, sort_idx] = sort(lfp_times);
            lfp_array = lfp_array(sort_idx);
        end
        
    catch ME
        warning('Error extracting LFP: %s', ME.message);
    end
end