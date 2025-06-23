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
    
    % find the first LFP row within tolerance
    dt = abs(lfpTimes - tdTime);
    idx = find(dt <= tol, 1);
    if isempty(idx)
        continue;
    end
    
    % split LFP channels, pick the one matching tdChan
    lfpChs = strsplit(BrainSenseLfp(idx).Channel, ',');
    if ~any(strcmp(lfpChs, tdChan))
        continue;
    end
    
    % grab the correct side's RateInHertz
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
        plot(stim_relative_time_sec, all_stim_amps, 'ro-', ...
             'LineWidth', 2, 'MarkerSize', 4, 'MarkerFaceColor', 'r');
        ylabel('Stimulation Amplitude (mA)', 'FontSize', 12, 'FontWeight', 'bold');
        ax = gca;
        ax.YColor = 'r';
        
        % Add legend
        legend({'Raw Signal', 'Stim Amplitude'}, 'Location', 'best');
    end
    
    % Add labels and title
    xlabel('Time (seconds)', 'FontSize', 12, 'FontWeight', 'bold');
    % Replace underscores with escaped underscores in the title
    escapedCombo = strrep(selectedCombo, '_', '\_');
    title(sprintf('Raw Time Domain Signal with Stimulation Amplitude: %s', escapedCombo), ...
          'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    % Store the current axes before any modifications
    current_ax = gca;
    
    % Detect and highlight dropped packet regions
    hold on;
    
    % Find gaps in the time series (more than 50ms between consecutive samples indicates dropped packets)
    time_diffs = diff(all_times);
    gap_threshold = milliseconds(50); % Normal sampling is 4ms, so 50ms indicates dropped packets
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
            patch([gap_start_time gap_end_time gap_end_time gap_start_time], ...
                  [ylims_left(1) ylims_left(1) ylims_left(2) ylims_left(2)], ...
                  [0.9 0.9 0.9], 'FaceAlpha', 0.5, 'EdgeColor', 'none', ...
                  'HandleVisibility', 'off');
            
            % Add text label for longer gaps (> 1 second)
            if gap_duration > 1
                text(mean([gap_start_time gap_end_time]), ylims_left(2)*0.9, ...
                     'Dropped Packets', 'HorizontalAlignment', 'center', ...
                     'FontSize', 10, 'Rotation', 90, 'Color', [0.5 0.5 0.5]);
            end
            
            fprintf('  Gap %d: %.1f - %.1f seconds (duration: %.2f s)\n', ...
                    i, gap_start_time, gap_end_time, gap_duration);
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

% Push updated table back to base workspace
assignin('base','BrainSenseTimeDomain',BrainSenseTimeDomain);
fprintf('\nProcessing complete!\n');