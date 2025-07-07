% 1) Prompt for JSON file and decode
[fileName, filePath] = uigetfile('*.json','Select a BrainSense JSON file');
if isequal(fileName,0)
    error('No file selected. Exiting.');
end
jsonText = fileread(fullfile(filePath,fileName));
data = jsondecode(jsonText);

% 2) Extract the two fields
if ~isfield(data,'BrainSenseTimeDomain') || ~isfield(data,'BrainSenseLfp')
    error('JSON missing required fields.');
end
rawTD = data.BrainSenseTimeDomain;
BrainSenseLfp = data.BrainSenseLfp;

% 3) Convert to table if needed
if ~istable(rawTD)
    BrainSenseTimeDomain = struct2table(rawTD);
else
    BrainSenseTimeDomain = rawTD;
end

% 4) Initialize StimRateHz
nTD = height(BrainSenseTimeDomain);
BrainSenseTimeDomain.StimRateHz = nan(nTD,1);

% 5) Parse timestamps once
tdTimes  = datetime(BrainSenseTimeDomain.FirstPacketDateTime, 'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
nLFP = numel(BrainSenseLfp);
lfpTimes = datetime({BrainSenseLfp.FirstPacketDateTime}, 'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');

% 6) Loop & match ±1.5 s to get StimRateHz
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

% 7) Group data by Channel and StimRateHz
fprintf('\nGrouping data by channel and stimulation rate...\n');

% Create unique combinations of Channel and StimRateHz
channels = BrainSenseTimeDomain.Channel;
stimRates = BrainSenseTimeDomain.StimRateHz;

% Create a combined identifier
channelRateCombos = {};
comboIndices = {};  % Store indices for each combo
for i = 1:nTD
    if ~isnan(stimRates(i))
        combo = sprintf('%s @ %.1f Hz', channels{i}, stimRates(i));
        
        % Check if this combo already exists
        existingIdx = find(strcmp(channelRateCombos, combo));
        if isempty(existingIdx)
            % New combo
            channelRateCombos{end+1} = combo;
            comboIndices{end+1} = i;
        else
            % Add to existing combo
            comboIndices{existingIdx} = [comboIndices{existingIdx}, i];
        end
    end
end

% Display available combinations
fprintf('\nAvailable Channel/StimRate combinations:\n');
for i = 1:length(channelRateCombos)
    fprintf('%d. %s (%d segments)\n', i, channelRateCombos{i}, length(comboIndices{i}));
end

% Prompt user for selection
selection = input('\nEnter the number of the combination you want to analyze: ');

if selection < 1 || selection > length(channelRateCombos)
    error('Invalid selection.');
end

selectedCombo = channelRateCombos{selection};
selectedIndices = comboIndices{selection};
fprintf('\nSelected: %s\n', selectedCombo);
fprintf('Processing %d segments...\n', length(selectedIndices));

% 8) Extract stimulation amplitude data for selected channel
fprintf('\nExtracting stimulation amplitude data...\n');

% Parse the selected channel to determine which side (LEFT/RIGHT)
if contains(selectedCombo, '_LEFT')
    targetSide = 'Left';
elseif contains(selectedCombo, '_RIGHT')
    targetSide = 'Right';
else
    error('Cannot determine side from channel name');
end

% Extract the base channel name (without LEFT/RIGHT)
channelBase = strrep(selectedCombo, ' @ ', '_');
channelBase = strrep(channelBase, sprintf('_%.1f Hz', stimRates(selectedIndices(1))), '');

fprintf('Looking for stim amp data for %s side...\n', targetSide);

% Collect all stimulation amplitude data for matching segments
all_stim_times = [];
all_stim_amps = [];

for idx = selectedIndices
    % Find matching LFP entry within tolerance
    tdTime = tdTimes(idx);
    dt = abs(lfpTimes - tdTime);
    lfpIdx = find(dt <= tol, 1);
    
    if isempty(lfpIdx)
        continue;
    end
    
    % Check if this LFP entry contains our channel
    lfpChannels = strsplit(BrainSenseLfp(lfpIdx).Channel, ',');
    if ~any(contains(lfpChannels, channelBase))
        continue;
    end
    
    % Extract LfpData
    lfpData = BrainSenseLfp(lfpIdx).LfpData;
    
    % Use the TimeDomain start time for alignment
    segmentStartTime = tdTimes(idx);
    
    % Add initial point at t=0 with stim amp = 0
    all_stim_times = [all_stim_times; segmentStartTime];
    all_stim_amps = [all_stim_amps; 0];
    
    % Process each element in the LfpData array
    for i = 1:length(lfpData)
        % Check if the target side exists
        if isfield(lfpData(i), targetSide)
            sideData = lfpData(i).(targetSide);
            
            % Extract mA value if it exists
            if isfield(sideData, 'mA')
                % Each LfpData point is 0.5 seconds apart (2 Hz sampling)
                % First measurement is at 0.5s, second at 1.0s, etc.
                stimTime = segmentStartTime + seconds(i * 0.5);
                stimAmp = sideData.mA;
                
                all_stim_times = [all_stim_times; stimTime];
                all_stim_amps = [all_stim_amps; stimAmp];
            end
        end
    end
end

% Sort stimulation data by time
if ~isempty(all_stim_times)
    [all_stim_times, sort_idx] = sort(all_stim_times);
    all_stim_amps = all_stim_amps(sort_idx);
    fprintf('Found %d stimulation amplitude data points\n', length(all_stim_amps));
else
    fprintf('No stimulation amplitude data found\n');
end

% 9) Plot raw time domain data and stimulation amplitude
figure('Position', [100, 100, 1400, 800]);

% Concatenate all time domain data for selected segments
all_times = [];
all_signals = [];

for idx = selectedIndices
    % Get the time domain data
    if ismember('TimeDomainData', BrainSenseTimeDomain.Properties.VariableNames) && ~isempty(BrainSenseTimeDomain.TimeDomainData{idx})
        
        signal = BrainSenseTimeDomain.TimeDomainData{idx};
        
        % Ensure signal is a column vector
        if size(signal, 2) > size(signal, 1)
            signal = signal';
        end
        
        % Get the start time for this segment
        startTime = tdTimes(idx);
        
        % Create time vector - each sample is 4ms apart
        nSamples = length(signal);
        sample_interval_ms = 4; % milliseconds
        
        % Generate timestamps for each sample
        sampleTimes = startTime + milliseconds((0:nSamples-1) * sample_interval_ms);
        
        % Append to our arrays
        all_times = [all_times; sampleTimes(:)];
        all_signals = [all_signals; signal(:)];
    end
end

% Sort by time to ensure proper ordering
[all_times, sort_idx] = sort(all_times);
all_signals = all_signals(sort_idx);

% Convert to relative time in seconds from first sample
if ~isempty(all_times)
    start_time = all_times(1);
    relative_time_sec = seconds(all_times - start_time);
    
    % Create dual y-axis plot
    yyaxis left
    
    % Plot the raw signal
    plot(relative_time_sec, all_signals, 'b-', 'LineWidth', 0.5);
    ylabel('Signal Amplitude', 'FontSize', 12, 'FontWeight', 'bold');
    ax = gca;
    ax.YColor = 'b';
    
    % Plot stimulation amplitude on right y-axis if available
    if ~isempty(all_stim_times)
        yyaxis right
        
        % Convert stim times to relative seconds
        stim_relative_time_sec = seconds(all_stim_times - start_time);
        
        % Plot stimulation amplitude
        plot(stim_relative_time_sec, all_stim_amps, 'ro-', 'LineWidth', 2, 'MarkerSize', 4, 'MarkerFaceColor', 'r');
        ylabel('Stimulation Amplitude (mA)', 'FontSize', 12, 'FontWeight', 'bold');
        ax = gca;
        ax.YColor = 'r';
        
        % Add legend
        legend({'Raw Signal', 'Stim Amplitude'}, 'Location', 'best');
    end
    
    % Add labels and title
    xlabel('Time (seconds) | Relative to inital value of firstpacketdatetime', 'FontSize', 12, 'FontWeight', 'bold');
    % Replace underscores with escaped underscores in the title
    escapedCombo = strrep(selectedCombo, '_', '\_');
    title(sprintf('Raw Time Domain Signal with Stimulation Amplitude @ %s', escapedCombo), 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    % Store the current axes before any modifications
    current_ax = gca;
    
    % Detect and highlight dropped packet regions
    hold on;
    
    % Find gaps in the time series (more than 50ms between consecutive samples indicates dropped packets)
    time_diffs = diff(all_times);
    gap_threshold = milliseconds(50); % Normal sampling is 4ms, so 50ms indicates dropped packets, can change this later down the line if needed
    gap_starts = find(time_diffs > gap_threshold);
    
    if ~isempty(gap_starts)
        fprintf('\nDetected %d dropped packet regions:\n', length(gap_starts));
        
        % Get current y-axis limits for shading (using left axis)
        yyaxis left
        ylims_left = ylim;
        yyaxis right
        ylims_right = ylim;
        
        % Use left axis limits for patch placement
        yyaxis left
        
        for i = 1:length(gap_starts)
            gap_start_time = seconds(all_times(gap_starts(i)) - start_time);
            gap_end_time = seconds(all_times(gap_starts(i) + 1) - start_time);
            gap_duration = gap_end_time - gap_start_time;
            
            % Shade the dropped packet region
            patch([gap_start_time gap_end_time gap_end_time gap_start_time], [ylims_left(1) ylims_left(1) ylims_left(2) ylims_left(2)], [0.9 0.9 0.9], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            
            % Add text label for longer gaps (> 1 second)
            if gap_duration > 1
                text(mean([gap_start_time gap_end_time]), ylims_left(2)*0.9, 'Dropped Packets', 'HorizontalAlignment', 'center', 'FontSize', 10, 'Rotation', 90, 'Color', [0.5 0.5 0.5]);
            end
            
            fprintf('Gap %d: %.1f - %.1f seconds (duration: %.2f s)\n', i, gap_start_time, gap_end_time, gap_duration);
        end
    end
    
    hold off;
    
    % Align zero lines of both y-axes
    if numel(current_ax.YAxis) == 2
        % Switch to left axis and get limits
        yyaxis left
        ylim_left = ylim;
        
        % Switch to right axis and get limits
        yyaxis right
        ylim_right = ylim;
        
        % Calculate ranges and zero positions
        range_left = ylim_left(2) - ylim_left(1);
        range_right = ylim_right(2) - ylim_right(1);
        
        % Find where zero is positioned in each axis (as fraction from bottom)
        if ylim_left(1) < 0 && ylim_left(2) > 0
            zero_pos_left = -ylim_left(1) / range_left;
        elseif ylim_left(1) >= 0
            zero_pos_left = 0; % Zero is at or below the bottom
        else
            zero_pos_left = 1; % Zero is at or above the top
        end
        
        if ylim_right(1) < 0 && ylim_right(2) > 0
            zero_pos_right = -ylim_right(1) / range_right;
        elseif ylim_right(1) >= 0
            zero_pos_right = 0; % Zero is at or below the bottom
        else
            zero_pos_right = 1; % Zero is at or above the top
        end
        
        % Determine target zero position (use the one that shows more data)
        target_zero_pos = max(zero_pos_left, zero_pos_right);
        
        % Adjust left axis
        if target_zero_pos > 0 && target_zero_pos < 1
            % Calculate new limits to position zero correctly
            new_range_left = max(ylim_left(2) / (1 - target_zero_pos), -ylim_left(1) / target_zero_pos);
            new_ylim_left = [-new_range_left * target_zero_pos, new_range_left * (1 - target_zero_pos)];
            
            % Ensure we don't cut off any data
            new_ylim_left(1) = min(new_ylim_left(1), ylim_left(1));
            new_ylim_left(2) = max(new_ylim_left(2), ylim_left(2));
            
            yyaxis left
            ylim(new_ylim_left);
        end
        
        % Adjust right axis
        if target_zero_pos > 0 && target_zero_pos < 1
            % Calculate new limits to position zero correctly
            new_range_right = max(ylim_right(2) / (1 - target_zero_pos), -ylim_right(1) / target_zero_pos);
            new_ylim_right = [-new_range_right * target_zero_pos, new_range_right * (1 - target_zero_pos)];
            
            % Ensure we don't cut off any data
            new_ylim_right(1) = min(new_ylim_right(1), ylim_right(1));
            new_ylim_right(2) = max(new_ylim_right(2), ylim_right(2));
            
            yyaxis right
            ylim(new_ylim_right);
        end
    end
    
    % Add statistics
    fprintf('\nSignal Statistics:\n');
    fprintf('Total duration: %.2f seconds\n', relative_time_sec(end));
    fprintf('Total samples: %d\n', length(all_signals));
    fprintf('Signal range: [%.2f, %.2f]\n', min(all_signals), max(all_signals));
    fprintf('Raw signal sampling: %d Hz (%.1f ms between samples)\n', 1000/sample_interval_ms, sample_interval_ms);
    
    if ~isempty(all_stim_amps)
        fprintf('\nStimulation Amplitude Statistics:\n');
        fprintf('Stim amp range: [%.2f, %.2f] mA\n', min(all_stim_amps), max(all_stim_amps));
        fprintf('Stim amp sampling: 2 Hz (0.5 s between samples)\n');
        fprintf('Number of stim measurements: %d\n', length(all_stim_amps));
    end
    
else
    fprintf('No time domain data found for the selected combination.\n');
end

% 10) Compute gamma power vs stimulation amplitude
% this is the new stuff so might need to come back to this and get older versiosn in github

% Get the stimulation rate for the selected channel
selected_stim_rate = stimRates(selectedIndices(1)); % Get stim rate from first selected index

% Calculate gamma band as 1/2 the stim rate (±1 Hz)
if ~isnan(selected_stim_rate)
    gamma_center = selected_stim_rate / 2;
    gamma_band = [gamma_center - 1, gamma_center + 1];
    fprintf('Selected stimulation rate: %.1f Hz\n', selected_stim_rate);
    fprintf('Calculated gamma band: [%.1f, %.1f] Hz (targeting %.1f Hz subharmonic)\n', gamma_band(1), gamma_band(2), gamma_center);
else
    error('No valid stimulation rate found for selected channel');
end

% Parameters
epoch_duration_sec = 10;  % 10-second epochs, also maybe experiment wtih this
fs = 250;  % Sampling frequency (Hz)
% gamma_band is now calculated above based on stim rate

% Welch's method parameters
welch_window_sec = 2; % 2-second window
welch_overlap = 0.5; % 50% overlap
welch_window_samples = welch_window_sec * fs;
welch_overlap_samples = floor(welch_window_samples * welch_overlap);

fprintf('\n\n=== Gamma Power vs Stimulation Amplitude Analysis ===\n');
fprintf('Epoch duration: %d seconds\n', epoch_duration_sec);
fprintf('Gamma band: %.1f-%.1f Hz (targeting %.1f Hz subharmonic response)\n', gamma_band(1), gamma_band(2), gamma_center);

% Check if we have data
if isempty(all_times) || isempty(all_signals)
    fprintf('No time domain data available for gamma analysis.\n');
    return;
end

% Initialize storage for results
epoch_results = [];

% Convert epoch duration to samples
epoch_samples = epoch_duration_sec * fs;

% Find the total duration of the recording
total_duration_sec = seconds(all_times(end) - all_times(1));
num_epochs = floor(total_duration_sec / epoch_duration_sec);

fprintf('Total recording duration: %.1f seconds\n', total_duration_sec);
fprintf('Number of potential epochs: %d\n', num_epochs);

% Process each epoch
valid_epochs = 0;
for epoch_idx = 1:num_epochs
    % Define epoch boundaries
    epoch_start_sec = (epoch_idx - 1) * epoch_duration_sec;
    epoch_end_sec = epoch_idx * epoch_duration_sec;
    
    % Convert to datetime boundaries
    epoch_start_time = all_times(1) + seconds(epoch_start_sec);
    epoch_end_time = all_times(1) + seconds(epoch_end_sec);
    
    % Find samples within this epoch
    epoch_mask = (all_times >= epoch_start_time) & (all_times < epoch_end_time);
    epoch_times = all_times(epoch_mask);
    epoch_signal = all_signals(epoch_mask);
    
    if isempty(epoch_times)
        continue;
    end
    
    % Check for dropped packets in this epoch
    % Calculate time differences between consecutive samples
    time_diffs_ms = milliseconds(diff(epoch_times));
    expected_interval_ms = 4;  % Expected 4ms between samples at 250Hz
    tolerance_ms = 1;  % Allow 1ms tolerance
    
    % Check if any gap exceeds expected interval
    has_dropped_packets = any(time_diffs_ms > (expected_interval_ms + tolerance_ms));
    
    % Also check if we have the expected number of samples
    expected_samples = epoch_duration_sec * fs;
    sample_count_ok = abs(length(epoch_signal) - expected_samples) <= (0.01 * expected_samples); % 1% tolerance
    
    if has_dropped_packets || ~sample_count_ok
        fprintf('Epoch %d: Skipped (dropped packets or incomplete data)\n', epoch_idx);
        continue;
    end
    
    % Valid epoch - compute gamma power
    try
        % Compute power spectral density using Welch's method
        [pxx, f] = pwelch(epoch_signal, hamming(welch_window_samples), welch_overlap_samples, [], fs);
        
        % Extract gamma band power
        gamma_idx = (f >= gamma_band(1)) & (f <= gamma_band(2));
        gamma_power = mean(pxx(gamma_idx));  % Average power in gamma band
        
        % Compute average stimulation amplitude for this epoch
        if ~isempty(all_stim_times) && ~isempty(all_stim_amps)
            % Find stim measurements within this epoch
            stim_mask = (all_stim_times >= epoch_start_time) & (all_stim_times < epoch_end_time);
            epoch_stim_amps = all_stim_amps(stim_mask);
            
            if ~isempty(epoch_stim_amps)
                avg_stim_amp = mean(epoch_stim_amps);
            else
                % No stim data in this epoch - interpolate from nearest values
                % Find nearest stim measurements before and after epoch
                before_idx = find(all_stim_times < epoch_start_time, 1, 'last');
                after_idx = find(all_stim_times > epoch_end_time, 1, 'first');
                
                if ~isempty(before_idx) && ~isempty(after_idx)
                    % Linear interpolation
                    t1 = seconds(all_stim_times(before_idx) - all_times(1));
                    t2 = seconds(all_stim_times(after_idx) - all_times(1));
                    epoch_mid_time = epoch_start_sec + epoch_duration_sec/2;
                    
                    avg_stim_amp = interp1([t1, t2], [all_stim_amps(before_idx), all_stim_amps(after_idx)], epoch_mid_time, 'linear');
                elseif ~isempty(before_idx)
                    avg_stim_amp = all_stim_amps(before_idx);
                elseif ~isempty(after_idx)
                    avg_stim_amp = all_stim_amps(after_idx);
                else
                    % No stim data available
                    continue;
                end
            end
        else
            % No stimulation data available
            avg_stim_amp = 0;
        end
        
        % Store results
        epoch_results(end+1, :) = [avg_stim_amp, gamma_power];
        valid_epochs = valid_epochs + 1;
        
        if mod(epoch_idx, 10) == 0
            fprintf('Processed epoch %d/%d\n', epoch_idx, num_epochs);
        end
        
    catch ME
        fprintf('Epoch %d: Error computing gamma power - %s\n', epoch_idx, ME.message);
        continue;
    end
end

fprintf('\nValid epochs processed: %d/%d\n', valid_epochs, num_epochs);

% Create scatter plot of gamma power vs stimulation amplitude
if ~isempty(epoch_results)
    figure('Position', [100, 100, 800, 600]);
    
    % Extract data
    stim_amps = epoch_results(:, 1);
    gamma_powers = epoch_results(:, 2);
    
    % Create scatter plot
    scatter(stim_amps, gamma_powers, 100, 'filled', 'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    
    % Add labels and title
    xlabel('Stimulation Amplitude (mA)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Gamma Power (μV²/Hz)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('Gamma Power (%d-%d Hz, 65 Hz Entrained) vs Stimulation Amplitude', gamma_band(1), gamma_band(2)), 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    % Add trend line if there are enough points
    if length(stim_amps) >= 3
        % Fit linear regression
        p = polyfit(stim_amps, gamma_powers, 1);
        
        % Generate trend line
        x_trend = linspace(min(stim_amps), max(stim_amps), 100);
        y_trend = polyval(p, x_trend);
        
        hold on;
        plot(x_trend, y_trend, 'r-', 'LineWidth', 2);
        
        % Calculate R-squared
        y_fit = polyval(p, stim_amps);
        SS_res = sum((gamma_powers - y_fit).^2);
        SS_tot = sum((gamma_powers - mean(gamma_powers)).^2);
        R_squared = 1 - (SS_res / SS_tot);
        
        % Add legend with statistics
        legend({'Data points', sprintf('Linear fit (R² = %.3f)', R_squared)}, 'Location', 'best', 'FontSize', 10);
        
        % Add text with equation
        equation_text = sprintf('y = %.2e·x + %.2e', p(1), p(2));
        text(0.05, 0.95, equation_text, 'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'white', 'EdgeColor', 'black');
    end
    
    hold off;
    
    % Print summary statistics
    fprintf('\n=== Results Summary ===\n');
    fprintf('Number of valid epochs: %d\n', size(epoch_results, 1));
    fprintf('Stimulation amplitude range: [%.2f, %.2f] mA\n', min(stim_amps), max(stim_amps));
    fprintf('Gamma power range: [%.2e, %.2e] μV²/Hz\n', min(gamma_powers), max(gamma_powers));
    
    if length(stim_amps) >= 3
        % Calculate correlation coefficient using corrcoef (built-in MATLAB)
        R_matrix = corrcoef(stim_amps, gamma_powers);
        correlation_coeff = R_matrix(1,2);
        fprintf('Linear correlation coefficient: %.3f\n', correlation_coeff);
        fprintf('Slope: %.2e μV²/Hz per mA\n', p(1));
    end
    
    % Create results table for export
    gamma_stim_table = table(stim_amps, gamma_powers, 'VariableNames', {'StimAmplitude_mA', 'GammaPower_uV2_Hz'});
    
    % Save to workspace
    assignin('base', 'gamma_stim_results', gamma_stim_table);
    assignin('base', 'epoch_results', epoch_results);
    
    fprintf('\nResults saved to workspace as ''gamma_stim_results'' (table) and ''epoch_results'' (array)\n');
    
else
    fprintf('\nNo valid epochs found for analysis.\n');
    fprintf('This could be due to:\n');
    fprintf('  - Too many dropped packets in the recording\n');
    fprintf('  - Recording duration shorter than epoch length\n');
    fprintf('  - Missing stimulation amplitude data\n');
end

% Making a heatmap if there are enough data points
if size(epoch_results, 1) >= 20
    figure('Position', [950, 100, 600, 500]);
    
    % Create 2D histogram
    num_bins = 10;
    stim_edges = linspace(min(stim_amps), max(stim_amps), num_bins);
    gamma_edges = linspace(min(gamma_powers), max(gamma_powers), num_bins);
    
    % Count occurrences in each bin
    [N, ~, ~] = histcounts2(stim_amps, gamma_powers, stim_edges, gamma_edges);
    
    % Create heatmap
    imagesc(stim_edges(1:end-1), gamma_edges(1:end-1), N');
    axis xy;
    colorbar;
    colormap(hot);
    
    xlabel('Stimulation Amplitude (mA)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Gamma Power (μV²/Hz)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Density Heatmap: Gamma Power vs Stimulation', 'FontSize', 14, 'FontWeight', 'bold');
    
    fprintf('\nDensity heatmap created (requires 20+ data points)\n');
end

fprintf('\n=== Analysis Complete ===\n');