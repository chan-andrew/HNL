%% 1) Prompt for JSON file and decode
[fileName, filePath] = uigetfile('*.json','Select a BrainSense JSON file');
if isequal(fileName,0)
    error('No file selected. Exiting.');
end
jsonText = fileread(fullfile(filePath,fileName));
data = jsondecode(jsonText);
%% 2) Extract the two fields
if ~isfield(data,'BrainSenseTimeDomain') || ~isfield(data,'BrainSenseLfp')
    error('JSON missing required fields.');
end
rawTD = data.BrainSenseTimeDomain;
BrainSenseLfp = data.BrainSenseLfp;
%% 3) Convert to table if needed
if ~istable(rawTD)
    BrainSenseTimeDomain = struct2table(rawTD);
else
    BrainSenseTimeDomain = rawTD;
end
%% 4) Initialize StimRateHz
nTD = height(BrainSenseTimeDomain);
BrainSenseTimeDomain.StimRateHz = nan(nTD,1);
%% 5) Parse timestamps once
tdTimes  = datetime(BrainSenseTimeDomain.FirstPacketDateTime, 'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
nLFP = numel(BrainSenseLfp);
lfpTimes = datetime({BrainSenseLfp.FirstPacketDateTime}, 'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
%% 6) Loop & match ±1.5 s
tol = seconds(1.5); % can change this down the line if the buffer needs to be tighter
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
%% 8) Build a "ThreePart" column of openable sub‐tables containing DateTime (placeholders for StimAmp & GammaPower)
nTD = height(BrainSenseTimeDomain);
% convert the original ISO strings to datetime once
DateTime = datetime(BrainSenseTimeDomain.FirstPacketDateTime, 'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
% preallocate cell array of subtables
threePartCells = cell(nTD,1);
for i = 1:nTD
    % each cell gets a 1×3 table:
    %   DateTime – actual timestamp
    %   StimAmp – NaN placeholder
    %   GammaPower – NaN placeholder
    threePartCells{i} = table( ...
        DateTime(i), ...
        nan, ...  % fill in with stimulation amplitude (mA) below
        nan, ...  % fill in with with gamma‐band power
        'VariableNames', {'DateTime','AvgStimAmp','GammaPower'} );
end
% attach it as a new column to your main table
BrainSenseTimeDomain.ThreePart = threePartCells;
%% 10) Push updated table back to base workspace
assignin('base','BrainSenseTimeDomain',BrainSenseTimeDomain);
%% 11) Debug version - Extract Stimulation Amplitude from BrainSenseLfp
fprintf('Starting stimulation amplitude extraction...\n');
fprintf('Number of TimeDomain rows: %d\n', nTD);
fprintf('Number of LFP entries: %d\n', nLFP);
for i = 1:nTD
    fprintf('\n--- Processing TimeDomain row %d ---\n', i);
    
    tdChan = BrainSenseTimeDomain.Channel{i};
    tdTime = tdTimes(i);
    fprintf('TD Channel: %s, TD Time: %s\n', tdChan, char(tdTime));
    
    % Find matching LFP entry within tolerance
    dt = abs(lfpTimes - tdTime);
    idx = find(dt <= tol, 1);
    
    if isempty(idx)
        fprintf('No matching LFP entry found within tolerance\n');
        continue;
    end
    
    fprintf('Found matching LFP entry at index %d\n', idx);
    fprintf('LFP Time: %s, Time difference: %.2f seconds\n', char(lfpTimes(idx)), seconds(dt(idx)));
    
    % Check if channels match
    lfpChannelStr = BrainSenseLfp(idx).Channel;
    lfpChs = strsplit(lfpChannelStr, ',');
    fprintf('LFP Channels: %s\n', lfpChannelStr);
    fprintf('Split LFP Channels: %s\n', strjoin(lfpChs, ', '));
    
    if ~any(strcmp(lfpChs, tdChan))
        fprintf('Channel %s not found in LFP channels\n', tdChan);
        continue;
    end
    
    fprintf('Channel match found!\n');
    
    % Determine which side we need
    if endsWith(tdChan, '_LEFT')
        targetSide = 'Left';
    elseif endsWith(tdChan, '_RIGHT')
        targetSide = 'Right';
    else
        fprintf('Unknown channel format: %s\n', tdChan);
        continue;
    end
    
    fprintf('Target side: %s\n', targetSide);
    
    % Check LfpData structure
    lfpData = BrainSenseLfp(idx).LfpData;
    fprintf('LfpData length: %d\n', length(lfpData));
    
    % Show structure of first LfpData entry
    if length(lfpData) > 0
        fprintf('LfpData(1) fields: %s\n', strjoin(fieldnames(lfpData(1)), ', '));
        if isfield(lfpData(1), targetSide)
            sideFields = fieldnames(lfpData(1).(targetSide));
            fprintf('LfpData(1).%s fields: %s\n', targetSide, strjoin(sideFields, ', '));
        else
            fprintf('Field %s not found in LfpData(1)\n', targetSide);
        end
    end
    
    % Extract mA values
    mA_values = [];
    for j = 1:length(lfpData)
        if isfield(lfpData(j), targetSide)
            sideData = lfpData(j).(targetSide);
            if isfield(sideData, 'mA')
                mA_val = sideData.mA;
                mA_values = [mA_values, mA_val];
                fprintf('LfpData(%d).%s.mA = %g\n', j, targetSide, mA_val);
            else
                fprintf('mA field not found in LfpData(%d).%s\n', j, targetSide);
            end
        else
            fprintf('Field %s not found in LfpData(%d)\n', targetSide, j);
        end
    end
    
    % Calculate average
    if ~isempty(mA_values)
        avgStimAmp = mean(mA_values);
        fprintf('Found %d mA values, average = %g\n', length(mA_values), avgStimAmp);
        
        % Update the ThreePart table
        BrainSenseTimeDomain.ThreePart{i}.AvgStimAmp = avgStimAmp;
        fprintf('Updated ThreePart table with AvgStimAmp = %g\n', avgStimAmp);
    else
        fprintf('No mA values found\n');
    end
end
%% 12) Robust Gamma Power Extraction (60-90 Hz) for all rows
fprintf('Starting gamma power extraction for all rows...\n');
% Define gamma frequency band
gamma_low = 60;   % Hz
gamma_high = 90;  % Hz
% Track statistics
successful_rows = 0;
failed_rows = 0;
for i = 1:nTD
    if mod(i, 10) == 0 || i <= 5  % Show progress every 10 rows, or first 5 rows
        fprintf('Processing row %d/%d\n', i, nTD);
    end
    
    try
        % Get the time domain signal data
        signal = [];
        
        % Check if TimeDomainData exists and has data
        if ismember('TimeDomainData', BrainSenseTimeDomain.Properties.VariableNames) && ...
           i <= length(BrainSenseTimeDomain.TimeDomainData) && ...
           ~isempty(BrainSenseTimeDomain.TimeDomainData{i})
            
            signal = BrainSenseTimeDomain.TimeDomainData{i};
        else
            if i <= 5  % Only show detailed messages for first few rows
                fprintf('No TimeDomainData for row %d\n', i);
            end
            failed_rows = failed_rows + 1;
            continue;
        end
        
        % Get sampling rate
        fs = BrainSenseTimeDomain.SampleRateInHz(i);
        
        % Validate sampling rate
        if isnan(fs) || fs <= 0 || fs > 10000  % reasonable upper bound
            if i <= 5
                fprintf('Invalid sampling rate for row %d: %g\n', i, fs);
            end
            failed_rows = failed_rows + 1;
            continue;
        end
        
        % Validate signal
        if ~isnumeric(signal) || isempty(signal)
            if i <= 5
                fprintf('Invalid signal for row %d\n', i);
            end
            failed_rows = failed_rows + 1;
            continue;
        end
        
        % Ensure signal is a column vector
        if size(signal, 2) > size(signal, 1)
            signal = signal';
        end
        
        % Remove NaN and Inf values
        signal = signal(isfinite(signal));
        
        % Check signal length after cleaning
        if length(signal) < 10  % Need at least 10 samples for meaningful analysis
            if i <= 5
                fprintf('Signal too short for row %d: %d samples\n', i, length(signal));
            end
            failed_rows = failed_rows + 1;
            continue;
        end
        
        % Check if gamma band is within Nyquist frequency
        nyquist_freq = fs / 2;
        if gamma_low >= nyquist_freq
            if i <= 5
                fprintf('Gamma band below Nyquist for row %d (fs=%g)\n', i, fs);
            end
            failed_rows = failed_rows + 1;
            continue;
        end
        
        % Adjust gamma_high if it exceeds Nyquist
        effective_gamma_high = min(gamma_high, nyquist_freq - 1);
        
        % Calculate power spectral density using pwelch
        [pxx, f] = pwelch(signal, [], [], [], fs);
        
        % Validate pwelch output
        if isempty(pxx) || isempty(f) || any(~isfinite(pxx)) || any(~isfinite(f))
            if i <= 5
                fprintf('Invalid pwelch output for row %d\n', i);
            end
            failed_rows = failed_rows + 1;
            continue;
        end
        
        % Find indices corresponding to gamma band
        gamma_indices = find(f >= gamma_low & f <= effective_gamma_high);
        
        if isempty(gamma_indices)
            if i <= 5
                fprintf('No gamma frequencies found for row %d\n', i);
            end
            failed_rows = failed_rows + 1;
            continue;
        end
        
        % Calculate gamma power as the mean power in the gamma band
        gamma_pxx = pxx(gamma_indices);
        gamma_power = mean(gamma_pxx);
        
        % Validate gamma power
        if ~isfinite(gamma_power) || gamma_power <= 0
            if i <= 5
                fprintf('Invalid gamma power for row %d: %g\n', i, gamma_power);
            end
            failed_rows = failed_rows + 1;
            continue;
        end
        
        % Convert to dB
        gamma_power_dB = 10 * log10(gamma_power);
        
        % Final validation of dB value
        if ~isfinite(gamma_power_dB)
            if i <= 5
                fprintf('Invalid gamma power dB for row %d: %g\n', i, gamma_power_dB);
            end
            failed_rows = failed_rows + 1;
            continue;
        end
        
        % Update the ThreePart table with gamma power
        BrainSenseTimeDomain.ThreePart{i}.GammaPower = gamma_power_dB;
        
        successful_rows = successful_rows + 1;
        
        if i <= 5  % Show details for first 5 rows
            fprintf('Row %d: SUCCESS - Gamma power = %.2f dB, Signal length = %d\n', ...
                i, gamma_power_dB, length(signal));
        end
        
    catch ME
        if i <= 5
            fprintf('Error processing row %d: %s\n', i, ME.message);
        end
        failed_rows = failed_rows + 1;
        continue;
    end
end
fprintf('\nGamma power extraction complete!\n');
fprintf('Successful rows: %d/%d\n', successful_rows, nTD);
fprintf('Failed rows: %d/%d\n', failed_rows, nTD);
%% 13) Push updated table back to base workspace
assignin('base','BrainSenseTimeDomain',BrainSenseTimeDomain);
fprintf('Updated BrainSenseTimeDomain saved to workspace.\n');
%% 14) Extract data from ThreePart tables and create plots
fprintf('\nExtracting data for plotting...\n');
% Initialize arrays to store extracted data
all_datetime = [];
all_stimamp = [];
all_gammapower = [];
% Extract data from each ThreePart subtable
for i = 1:nTD
    if ~isempty(BrainSenseTimeDomain.ThreePart{i})
        subTable = BrainSenseTimeDomain.ThreePart{i};
        
        % Extract datetime
        dt = subTable.DateTime;
        all_datetime = [all_datetime; dt];
        
        % Extract stim amp (handle NaN values)
        stimamp = subTable.AvgStimAmp;
        all_stimamp = [all_stimamp; stimamp];
        
        % Extract gamma power (handle NaN values)
        gammapower = subTable.GammaPower;
        all_gammapower = [all_gammapower; gammapower];
    end
end
% Convert to relative time in minutes
if ~isempty(all_datetime)
    % Sort by datetime to ensure proper ordering
    [all_datetime, sort_idx] = sort(all_datetime);
    all_stimamp = all_stimamp(sort_idx);
    all_gammapower = all_gammapower(sort_idx);
    
    % Calculate relative time in minutes from first timestamp
    start_time = all_datetime(1);
    relative_time_min = minutes(all_datetime - start_time);
    
    % Create figure with dual y-axes
    figure('Position', [100, 100, 1200, 600]);
    
    % Create left y-axis for stimulation amplitude
    yyaxis left
    
    % Plot stimulation amplitude (filter out NaN values for plotting)
    valid_stim_idx = ~isnan(all_stimamp);
    if any(valid_stim_idx)
        plot(relative_time_min(valid_stim_idx), all_stimamp(valid_stim_idx), ...
             'o-', 'LineWidth', 2, 'MarkerSize', 6, 'Color', [0.0, 0.4, 0.8]);
        ylabel('Stimulation Amplitude (mA)', 'FontSize', 12, 'FontWeight', 'bold');
    else
        ylabel('Stimulation Amplitude (mA)', 'FontSize', 12, 'FontWeight', 'bold');
    end
    
    % Set left y-axis color
    ax = gca;
    ax.YColor = [0.0, 0.4, 0.8];
    
    % Create right y-axis for gamma power
    yyaxis right
    
    % Plot gamma power (filter out NaN values for plotting)
    valid_gamma_idx = ~isnan(all_gammapower);
    if any(valid_gamma_idx)
        plot(relative_time_min(valid_gamma_idx), all_gammapower(valid_gamma_idx), ...
             's-', 'LineWidth', 2, 'MarkerSize', 6, 'Color', [0.8, 0.2, 0.2]);
        ylabel('Gamma Power (dB)', 'FontSize', 12, 'FontWeight', 'bold');
    else
        ylabel('Gamma Power (dB)', 'FontSize', 12, 'FontWeight', 'bold');
    end
    
    % Set right y-axis color
    ax = gca;
    ax.YColor = [0.8, 0.2, 0.2];
    
    % Set x-axis properties
    xlabel('Relative Time (minutes)', 'FontSize', 12, 'FontWeight', 'bold');
    title('BrainSense Data: Stimulation Amplitude and Gamma Power Over Time', ...
          'FontSize', 14, 'FontWeight', 'bold');
    
    % Add grid
    grid on;
    
    % Add legend
    legend_entries = {};
    if any(valid_stim_idx)
        legend_entries{end+1} = 'Stimulation Amplitude (mA)';
    end
    if any(valid_gamma_idx)
        legend_entries{end+1} = 'Gamma Power (dB)';
    end
    
    if ~isempty(legend_entries)
        legend(legend_entries, 'Location', 'best', 'FontSize', 10);
    end
    
    %% 14a) Align Y-Axis Zero Lines [NEW CODE]
    fprintf('\nAligning Y-axis zero lines...\n');
    
    ax = gca; % Get current axes
    
    % Proceed only if there are two y-axes
    if numel(ax.YAxis) == 2
        ylim_left_orig = ax.YAxis(1).Limits;
        ylim_right_orig = ax.YAxis(2).Limits;
        
        % --- Calculate target proportions ---
        % Proportion of the axis required for the part > 0
        range_L = ylim_left_orig(2) - ylim_left_orig(1);
        prop_L_pos = (range_L > 0) * ylim_left_orig(2) / range_L;
        
        range_R = ylim_right_orig(2) - ylim_right_orig(1);
        prop_R_pos = (range_R > 0) * ylim_right_orig(2) / range_R;
        
        % Clamp proportions to handle data that is entirely positive or negative
        prop_L_pos = max(0, min(1, prop_L_pos));
        prop_R_pos = max(0, min(1, prop_R_pos));
        
        % The target proportion for the positive side is the max of the two
        target_prop_pos = max(prop_L_pos, prop_R_pos);
        % The target proportion for the negative side is also the max of the two
        target_prop_neg = max(1 - prop_L_pos, 1 - prop_R_pos);
        
        % Rescale proportions so their sum is 1
        total_prop = target_prop_pos + target_prop_neg;
        if total_prop > 0
            final_prop_neg = target_prop_neg / total_prop;
        else
            final_prop_neg = 0.5; % Default for empty/flat data
        end
        
        % --- Adjust Left Y-Axis ---
        if final_prop_neg > 0 && final_prop_neg < 1
            % Calculate the required total range based on the most extreme data point
            range_from_pos = ylim_left_orig(2) / (1 - final_prop_neg);
            range_from_neg = -ylim_left_orig(1) / final_prop_neg;
            new_range_L = max(range_from_pos, range_from_neg);
            
            % Set new limits based on the calculated range and proportion
            ax.YAxis(1).Limits = [-new_range_L * final_prop_neg, new_range_L * (1 - final_prop_neg)];
        end
        
        % --- Adjust Right Y-Axis ---
        if final_prop_neg > 0 && final_prop_neg < 1
            % Calculate the required total range for the right axis
            range_from_pos = ylim_right_orig(2) / (1 - final_prop_neg);
            range_from_neg = -ylim_right_orig(1) / final_prop_neg;
            new_range_R = max(range_from_pos, range_from_neg);
            
            % Set new limits for the right axis
            ax.YAxis(2).Limits = [-new_range_R * final_prop_neg, new_range_R * (1 - final_prop_neg)];
        end
        
        fprintf('Zero lines aligned.\n');
    else
        fprintf('Skipping zero line alignment (only one Y-axis present).\n');
    end
    
    % Print summary statistics
    fprintf('\nPlotting Summary:\n');
    fprintf('Total data points: %d\n', length(all_datetime));
    fprintf('Time range: %.2f minutes\n', max(relative_time_min));
    
    if any(valid_stim_idx)
        fprintf('Valid stimulation amplitude points: %d\n', sum(valid_stim_idx));
        fprintf('Stim amp range: %.2f - %.2f mA\n', min(all_stimamp(valid_stim_idx)), max(all_stimamp(valid_stim_idx)));
    else
        fprintf('No valid stimulation amplitude data found\n');
    end
    
    if any(valid_gamma_idx)
        fprintf('Valid gamma power points: %d\n', sum(valid_gamma_idx));
        fprintf('Gamma power range: %.2f - %.2f dB\n', min(all_gammapower(valid_gamma_idx)), max(all_gammapower(valid_gamma_idx)));
    else
        fprintf('No valid gamma power data found\n');
    end
    
else
    fprintf('No datetime data found for plotting\n');
end

%% 15) Create Scatter Plot of Stimulation Amplitude vs. Gamma Power [NEW]

fprintf('\nCreating relationship plot: Stim Amp vs. Gamma Power...\n');

% Find the indices where BOTH variables have valid data (are not NaN)
valid_indices = ~isnan(all_stimamp) & ~isnan(all_gammapower);

% Extract only the valid data pairs
stim_amp_data = all_stimamp(valid_indices);
gamma_power_data = all_gammapower(valid_indices);

% Check if there is any data to plot
if isempty(stim_amp_data)
    fprintf('No overlapping data points found to create a relationship plot.\n');
else
    % Create a new figure for the scatter plot
    figure('Name', 'Stimulation vs. Gamma Power Relationship');

    % Create the scatter plot
    scatter(stim_amp_data, gamma_power_data, 50, 'filled', 'MarkerFaceColor', [0.2, 0.2, 0.7]);
    
    % --- Optional: Add a line of best fit (linear regression) ---
    hold on;
    % Calculate the coefficients for a 1st degree polynomial (a line)
    p = polyfit(stim_amp_data, gamma_power_data, 1); 
    % Evaluate the polynomial at the x-data points to get the y-values for the line
    y_fit = polyval(p, stim_amp_data);
    % Plot the line of best fit
    plot(stim_amp_data, y_fit, 'r-', 'LineWidth', 2);
    hold off;
    % --- End of optional section ---

    % Add labels, title, and grid for clarity
    xlabel('Stimulation Amplitude (mA)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Gamma Power (dB)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Relationship Between Gamma Power and Stimulation Amplitude', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    % Add a legend if you included the line of best fit
    legend('Data Points', 'Linear Fit', 'Location', 'best');
    
    fprintf('Relationship plot created with %d data points.\n', numel(stim_amp_data));
end
fprintf('\nProcessing complete!\n');