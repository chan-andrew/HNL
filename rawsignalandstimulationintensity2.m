function final = rawsignalandstimulationintensity2()

    % 1) Prompt for JSON file
    [fileName, filePath] = uigetfile('*.json','Select BrainSense JSON file');
    if isequal(fileName,0)
        fprintf('No file selected. Exiting.\n');
        final = [];
        return;
    end
    fullFile = fullfile(filePath,fileName);

    % 2) Load and decode JSON
    try
        jsonText = fileread(fullFile);
        data = jsondecode(jsonText);
    catch ME
        error('Error reading or decoding JSON: %s', ME.message);
    end

    % 3) Build timestamps for TD and LFP
    maxSeq = 255;
    TD = buildTimeStamps(data.BrainSenseTimeDomain, maxSeq, 'TD');
    LFP = buildTimeStamps(data.BrainSenseLfp, maxSeq, 'LFP');

    acc_orig_Lstim_times_utc = [];
    acc_orig_Lstim_values = [];
    acc_orig_Rstim_times_utc = [];
    acc_orig_Rstim_values = [];

    % 4) Align channels by FirstPacketDateTime

    tdFP = datetime({TD.FirstPacketDateTime},  'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
    lfpFP = datetime({LFP.FirstPacketDateTime}, 'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
    idxMap = zeros(size(lfpFP));
    for i = 1:numel(lfpFP)
        if ~isempty(tdFP) && ~isempty(lfpFP) && ~isempty(lfpFP(i)) && ~isempty(tdFP) % Ensure not empty before calculation
             if (isdatetime(lfpFP(i)) && ~isnat(lfpFP(i))) % Check if lfpFP(i) is a valid datetime
                [~,idxMap(i)] = min(abs(tdFP - lfpFP(i)));
            else
                idxMap(i) = NaN; % Handle if lfpFP(i) is not a valid datetime
            end
        else
            idxMap(i) = NaN; % Or handle error appropriately
        end
    end
    if numel(unique(idxMap(~isnan(idxMap)))) ~= numel(idxMap(~isnan(idxMap)))
        warning('One or more LFP channels align to the same TD channel based on FirstPacketDateTime proximity.');
    end

    % 5) Extract and assemble final table
    final = table();
    for j = 1:numel(LFP) % Loop through each LFP stream
        % Channel pairing from LFP(j).Channel string
        chanPair = split(LFP(j).Channel, ',');
        chanLeft = ''; chanRight = ''; % Initialize
        if numel(chanPair) >= 1, chanLeft  = strtrim(chanPair{1}); end
        if numel(chanPair) >= 2, chanRight = strtrim(chanPair{2}); end

        % Packet date string from LFP(j) to match with TD streams
        fpStr = LFP(j).FirstPacketDateTime;

        % Find matching TD stream indices based on Channel name and FirstPacketDateTime
        idxL = []; idxR = [];
        if ~isempty(chanLeft), idxL = find(strcmp({TD.Channel}, chanLeft) & strcmp({TD.FirstPacketDateTime}, fpStr)); end
        if ~isempty(chanRight), idxR = find(strcmp({TD.Channel}, chanRight) & strcmp({TD.FirstPacketDateTime}, fpStr)); end

        % Check if current LFP(j) stream contains stimulation mA data in its LfpData
        hasLeftStimInLFPData = false; hasRightStimInLFPData = false;
        current_maL_original = []; current_maR_original = []; % mA values from this LFP(j)

        if isfield(LFP(j), 'LfpData') && ~isempty(LFP(j).LfpData)
            % Check for Left stimulation data
            if isfield(LFP(j).LfpData(1), 'Left') && isfield(LFP(j).LfpData(1).Left, 'mA')
                hasLeftStimInLFPData = true;
                current_maL_original = arrayfun(@(pkt) pkt.Left.mA, LFP(j).LfpData, 'UniformOutput', true);
            end
            % Check for Right stimulation data
            if isfield(LFP(j).LfpData(1), 'Right') && isfield(LFP(j).LfpData(1).Right, 'mA')
                hasRightStimInLFPData = true;
                current_maR_original = arrayfun(@(pkt) pkt.Right.mA, LFP(j).LfpData, 'UniformOutput', true);
            end

            % If we found mA data, get their original timestamps
            if (hasLeftStimInLFPData || hasRightStimInLFPData) && isfield(LFP(j).LfpData(1), 'TicksInMs')
                packet_ticks = arrayfun(@(pkt) pkt.TicksInMs, LFP(j).LfpData);
                t0_LFP_stream = datetime(LFP(j).FirstPacketDateTime, 'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');

                if ~isempty(packet_ticks)
                    % Timestamps for each packet/event in LFP(j).LfpData
                    original_event_times_for_LFPj = t0_LFP_stream + milliseconds(packet_ticks - packet_ticks(1));

                    if hasLeftStimInLFPData && ~isempty(current_maL_original)
                        acc_orig_Lstim_times_utc = [acc_orig_Lstim_times_utc; original_event_times_for_LFPj(:)];
                        acc_orig_Lstim_values = [acc_orig_Lstim_values; current_maL_original(:)];
                    end
                    if hasRightStimInLFPData && ~isempty(current_maR_original)
                        acc_orig_Rstim_times_utc = [acc_orig_Rstim_times_utc; original_event_times_for_LFPj(:)];
                        acc_orig_Rstim_values = [acc_orig_Rstim_values; current_maR_original(:)];
                    end
                end
            end
        end

        if ~(numel(idxL)==1 && numel(idxR)==1)
            continue; % Skip to next LFP stream if TD channels not matched for table construction
        end

        % If here, idxL and idxR are valid, corresponding to TD streams
        tvec = TD(idxL).TimeStamps; % High-rate timestamps from the matched TD stream
        rawLeft = TD(idxL).TimeDomainData(:);
        rawRight = TD(idxR).TimeDomainData(:); % Assumes idxR is also valid and corresponds to rawRight

        % Stimulation currents for the table (from current_maL_original, current_maR_original if LFP(j) had them)
        % These are the low-rate (e.g., 2Hz) stimulation values
        maL_for_table = current_maL_original;
        maR_for_table = current_maR_original;

        % Prepare stimL and stimR for the table using repelem (expand low-rate stim to high-rate tvec)
        stimL = nan(size(rawLeft));  % Default to NaN
        stimR = nan(size(rawRight)); % Default to NaN

        if ~isempty(maL_for_table)
            nPackets_L = numel(maL_for_table);
            Nsig_L = numel(rawLeft);
            if nPackets_L > 0 && Nsig_L > 0
                reps_L = ceil(Nsig_L / nPackets_L);
                stimL_full = repelem(maL_for_table, reps_L);
                stimL = stimL_full(1:Nsig_L); % Trim to exact length of raw data
            end
        end

        if ~isempty(maR_for_table)
            nPackets_R = numel(maR_for_table);
            Nsig_R = numel(rawRight);
            if nPackets_R > 0 && Nsig_R > 0
                reps_R = ceil(Nsig_R / nPackets_R);
                stimR_full = repelem(maR_for_table, reps_R);
                stimR = stimR_full(1:Nsig_R); % Trim to exact length of raw data
            end
        end
        
        % Build table row for this channel pair and append to final table
        T = table(tvec, rawLeft, stimL, rawRight, stimR, ...
            'VariableNames', {'Time','LeftRaw','LeftStim_mA','RightRaw','RightStim_mA'});
        final = [final; T];
    end % End of for j = 1:numel(LFP)

    % 6) Assign to base workspace and return
    assignin('base','final', final);
    if ~isempty(final)
        disp('''final'' table created in workspace with 5 columns.');
        if any(diff(final.Time) < 0)
            disp('Sorting final table by time.');
            final = sortrows(final, 'Time');
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- BEGIN NEW SECTION: Analyze Time Differences for Dropped Packets/Gaps ---
        if height(final) > 1 % Ensure there's data and more than one point for diff
            fprintf('\nAnalyzing time differences in final.Time for potential gaps...\n');
        
            time_durations = diff(final.Time);       % Calculate duration objects
            time_diff_in_seconds = seconds(time_durations); % Convert to numeric seconds
        
            if ~isempty(time_diff_in_seconds)
                fprintf('Statistics for time_diff_in_seconds (gaps between consecutive samples):\n');
                fprintf('Min diff: %.6f s\n', min(time_diff_in_seconds));
                fprintf('Max diff: %.6f s\n', max(time_diff_in_seconds));
                fprintf('Mean diff: %.6f s\n', mean(time_diff_in_seconds));
                fprintf('Median diff: %.6f s (This is often your nominal sampling interval)\n', median(time_diff_in_seconds));
        
                % Estimate nominal sampling interval from the median difference
                nominal_interval = median(time_diff_in_seconds);
                if nominal_interval > 0
                    fprintf('Estimated nominal sampling frequency: %.2f Hz\n', 1/nominal_interval);
                end
        
                % Define a threshold for what constitutes a "significant" gap
                % Example: A gap is "significant" if it's more than 1.5 times the nominal interval. 
                % Adjust this multiplier or use an absolute threshold (e.g., 0.1 seconds) as needed.
                if nominal_interval > 0
                    gap_threshold_multiplier = 1.5;
                    gap_threshold = nominal_interval * gap_threshold_multiplier;
                else % Fallback if nominal_interval is zero or negative (should not happen with sorted positive time)
                    gap_threshold = 0.01; % Default to a small absolute threshold like 10ms
                end
        
                significant_gaps_indices = find(time_diff_in_seconds > gap_threshold);
        
                if ~isempty(significant_gaps_indices)
                    fprintf('\nFound %d significant time gaps (exceeding %.6f s):\n', numel(significant_gaps_indices), gap_threshold);
                    % Display details for the first few gaps (e.g., up to 10)
                    for i = 1:min(numel(significant_gaps_indices), 10)
                        idx = significant_gaps_indices(i);
                        fprintf('  Gap of %.6f s after final.Time(%d) (between original timestamps %s and %s)\n', ...
                            time_diff_in_seconds(idx), idx, ...
                            datestr(final.Time(idx),'yyyy-mm-dd HH:MM:SS.FFF'), ...
                            datestr(final.Time(idx+1),'yyyy-mm-dd HH:MM:SS.FFF'));
                    end
                    if numel(significant_gaps_indices) > 10
                        fprintf('... and %d more gaps.\n', numel(significant_gaps_indices)-10);
                    end
                else
                    fprintf('No significant time gaps detected based on the threshold of %.6f s.\n', gap_threshold);
                end
            else
                fprintf('No time differences to analyze (final.Time has less than 2 data points).\n');
            end
        
            % Plot the time differences
            figure('Name', 'Time Differences Between Consecutive Samples in final.Time');
            plot(time_diff_in_seconds);
            xlabel('Index of Difference (corresponds to final.Time index k and k+1)');
            ylabel('Time Difference (seconds)');
            title('Actual Time Difference Between Consecutive Samples in final.Time');
            grid on;
            if ~isempty(time_diff_in_seconds) && nominal_interval > 0
                % Adjust y-limits for better visualization, focusing on values around the nominal interval
                % and ensuring large gaps are visible.
                ylim_upper = max(nominal_interval * 10, max(time_diff_in_seconds) * 1.1);
                if ylim_upper == 0 || isnan(ylim_upper) % Handle edge cases
                    ylim_upper = 1;
                end
                if min(time_diff_in_seconds) < 0 % Should not happen if time is sorted
                    ylim_lower = min(time_diff_in_seconds) * 1.1;
                else
                    ylim_lower = 0;
                end
                ylim([ylim_lower, ylim_upper]);
            end
            fprintf('A new figure "Time Differences Between Consecutive Samples in final.Time" has been generated.\n\n');
        else
            fprintf('Skipping time difference analysis as ''final'' table has insufficient data (<= 1 row).\n\n');
        end
        % --- END NEW SECTION ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    else
        disp('''final'' table is empty. Check LFP and TD data and matching logic.');
    end

    % 7) Plotting Section
    fprintf('Plotting stimulation data...\n');
    h_fig = figure('Name', 'Stimulation Data Verification', 'NumberTitle', 'off'); % Get figure handle

    % Plot Left Stimulation
    h_ax1 = subplot(2,1,1); % Get axes handle for left plot
    hold(h_ax1, 'on');
    if ~isempty(acc_orig_Lstim_times_utc) && ~isempty(acc_orig_Lstim_values)
        % Plot unique original events to avoid overplotting if multiple LFP streams had same data
        [unique_orig_L_times, ia, ~] = unique(acc_orig_Lstim_times_utc);
        unique_orig_L_values = acc_orig_Lstim_values(ia);
        plot(h_ax1, unique_orig_L_times, unique_orig_L_values, 'bo-', 'DisplayName', 'Original Left Stim Events (LFP packets)', 'MarkerSize', 6);
        fprintf('Plotting %d unique original Left Stim event points.\n', numel(unique_orig_L_values));
    else
        fprintf('No original Left Stimulation event data collected to plot.\n');
    end

    if ~isempty(final) && any(strcmp(final.Properties.VariableNames, 'LeftStim_mA')) && ~all(isnan(final.LeftStim_mA))
        plot(h_ax1, final.Time, final.LeftStim_mA, 'r.-', 'DisplayName', 'Processed Left Stim (final table)', 'MarkerSize', 4);
        fprintf('Plotting processed Left Stim data from final table.\n');
    else
        fprintf('No valid Left Stimulation data in final table to plot.\n');
    end

    hold(h_ax1, 'off');
    title(h_ax1, 'Left Stimulation Intensity');
    xlabel(h_ax1, 'Time (UTC)');
    ylabel(h_ax1, 'Stimulation (mA)');
    legend(h_ax1, 'show');
    grid(h_ax1, 'on');
    % Initial datetick for full view (optional, as it will be applied after zoom)
    % datetick(h_ax1, 'x', 'HH:MM:SS.FFF', 'keeplimits'); 

    % Plot Right Stimulation
    h_ax2 = subplot(2,1,2); % Get axes handle for right plot
    hold(h_ax2, 'on');
    if ~isempty(acc_orig_Rstim_times_utc) && ~isempty(acc_orig_Rstim_values)
        [unique_orig_R_times, ia, ~] = unique(acc_orig_Rstim_times_utc);
        unique_orig_R_values = acc_orig_Rstim_values(ia);
        plot(h_ax2, unique_orig_R_times, unique_orig_R_values, 'go-', 'DisplayName', 'Original Right Stim Events (LFP packets)', 'MarkerSize', 6);
        fprintf('Plotting %d unique original Right Stim event points.\n', numel(unique_orig_R_values));
    else
        fprintf('No original Right Stimulation event data collected to plot.\n');
    end

    if ~isempty(final) && any(strcmp(final.Properties.VariableNames, 'RightStim_mA')) && ~all(isnan(final.RightStim_mA))
        plot(h_ax2, final.Time, final.RightStim_mA, 'm.-', 'DisplayName', 'Processed Right Stim (final table)', 'MarkerSize', 4);
        fprintf('Plotting processed Right Stim data from final table.\n');
    else
        fprintf('No valid Right Stimulation data in final table to plot.\n');
    end

    hold(h_ax2, 'off');
    title(h_ax2, 'Right Stimulation Intensity');
    xlabel(h_ax2, 'Time (UTC)');
    ylabel(h_ax2, 'Stimulation (mA)');
    legend(h_ax2, 'show');
    grid(h_ax2, 'on');
    % Initial datetick for full view (optional)
    % datetick(h_ax2, 'x', 'HH:MM:SS.FFF', 'keeplimits');

    sgtitle(h_fig, 'Stimulation Data Verification'); % Overall title for the figure

    % Zoom and Finalize Ticks
    if ~isempty(final) && ~isempty(final.Time) && isdatetime(final.Time) % Ensure final.Time is valid datetime
        % Focusing on this time stamp right now for testing reasons: '11-Mar-2025 14:40:30'
        if isempty(final.Time(1).TimeZone)
             warning('final.Time has no timezone. Assuming UTC for zooming.');
             timeOfInterest = datetime('11-Mar-2025 14:40:30', 'TimeZone', 'UTC');
        else
            timeOfInterest = datetime('11-Mar-2025 14:40:30', 'TimeZone', final.Time(1).TimeZone);
        end
        
        timeWindow = minutes(1); % Show a +/- 1 minute window

        % Apply zoom and datetick to Left Plot (h_ax1)
        xlim(h_ax1, [timeOfInterest - timeWindow, timeOfInterest + timeWindow]);
        datetick(h_ax1, 'x', 'HH:MM:SS.FFF', 'keeplimits'); % Apply datetick AFTER xlim
        
        % Apply zoom and datetick to Right Plot (h_ax2)
        xlim(h_ax2, [timeOfInterest - timeWindow, timeOfInterest + timeWindow]);
        datetick(h_ax2, 'x', 'HH:MM:SS.FFF', 'keeplimits'); % Apply datetick AFTER xlim

        fprintf('Zoomed plots to a window around 11-Mar-2025 14:40:30 (TZ: %s).\n', char(timeOfInterest.TimeZone));
    else
         fprintf('Final table is empty or has no valid Time data, skipping plot zoom. Applying default time ticks to full range if possible.\n');
         % Apply datetick to full range if zoom is skipped but plots were made
         if exist('h_ax1', 'var'), datetick(h_ax1, 'x', 'HH:MM:SS.FFF', 'keeplimits'); end
         if exist('h_ax2', 'var'), datetick(h_ax2, 'x', 'HH:MM:SS.FFF', 'keeplimits'); end
    end

    fprintf('Plotting complete.\n');
end % End of main function rawsignalandstimulationintensity

% Helper function: buildTimeStamps (your original function)
function S = buildTimeStamps(S, maxSeq, type)
    % buildTimeStamps builds S.TimeStamps for either 'TD' or 'LFP'
    for k = 1:numel(S)
        switch type
            case 'TD'
                ticks = parseNumeric(S(k).TicksInMses); % Note: field name TicksInMses
                sizes = parseNumeric(S(k).GlobalPacketSizes);
                if isfield(S(k),'GlobalSequences')
                    seqs = parseNumeric(S(k).GlobalSequences);
                else
                    seqs = 1:numel(ticks); % Simplified sequence if missing
                end
            case 'LFP'
                % TicksInMs is per packet in LfpData array
                if isfield(S(k),'LfpData') && ~isempty(S(k).LfpData) && isfield(S(k).LfpData(1),'TicksInMs')
                    rawTicks = arrayfun(@(pkt) pkt.TicksInMs, S(k).LfpData);
                    ticks = rawTicks(:);
                else
                    fprintf('Warning: LFP channel %d (%s) missing LfpData or TicksInMs. Timestamps may be incorrect.\n', k, S(k).Channel);
                    S(k).TimeStamps = []; % Cannot compute timestamps
                    continue; % Skip to next LFP stream
                end

                % GlobalPacketSizes might not be relevant if LFP data is 1 sample per LfpData entry
                if isfield(S(k),'GlobalPacketSizes')
                    sizes = parseNumeric(S(k).GlobalPacketSizes);
                    if numel(sizes) ~= numel(ticks) % Packet sizes should match number of LfpData packets
                        sizes = ones(size(ticks)); % Safer default if mismatch
                    end
                else
                    sizes = ones(size(ticks)); % Assume 1 sample per LfpData packet if GlobalPacketSizes is missing
                end

                if isfield(S(k),'GlobalSequences') % For LFP, this might not always be present or relevant
                    seqs = parseNumeric(S(k).GlobalSequences);
                elseif isfield(S(k),'LfpData') && ~isempty(S(k).LfpData) && isfield(S(k).LfpData(1),'Seq')
                    seqs = arrayfun(@(pkt) pkt.Seq, S(k).LfpData); % Sequence per packet
                else
                    seqs = 1:numel(ticks); % Simplified sequence
                end
            otherwise
                error('Unknown type ''%s''', type);
        end

        if isempty(ticks) % If no ticks, cannot proceed
            S(k).TimeStamps = [];
            continue;
        end

        fs   = S(k).SampleRateInHz;
        if fs == 0 % Avoid division by zero if SampleRateInHz is 0
            warning('SampleRateInHz is 0 for %s channel %d (%s). Timestamps may be incorrect.', type, k, S(k).Channel);
            dt = seconds(0); % Or handle as error / assign NaN
        else
            dt   = seconds(1/fs);
        end
        
        t0   = datetime(S(k).FirstPacketDateTime, ...
                     'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
        
        % tStart is the start time of each packet/segment
        tStart = t0 + milliseconds(ticks - ticks(1)); % Relative to the first tick in this stream

        % Check for dropped packets (optional, can be verbose)
        dseq = diff(seqs);
        dseq(dseq < 0) = dseq(dseq < 0) + maxSeq + 1; % Handle sequence wrap-around
        if any(dseq > 1) % Check for gaps greater than 1
            % warning('%s channel %d (%s): Potential dropped packets detected based on sequence numbers.', type, k, S(k).Channel);
        end

        segs = cell(numel(ticks),1); % One cell per packet/segment
        for i = 1:numel(ticks)
            if sizes(i) > 0
                segs{i} = (tStart(i) : dt : tStart(i) + dt*(sizes(i)-1)).';
            else
                segs{i} = []; % No samples if size is 0
            end
        end
        S(k).TimeStamps = vertcat(segs{:}); % Concatenate all detailed timestamps
    end
end

% Helper function: parseNumeric (your original function)
function v = parseNumeric(raw)
    % parseNumeric converts numeric, cell, or comma-delimited string to numeric column
    if isnumeric(raw)
        v = raw(:);
    elseif iscell(raw)
        try % Add try-catch for robustness if cell contains non-convertible strings
            v_temp = str2double(raw(:)); % Attempt conversion
            v = v_temp; % Assign if successful
        catch
            v = nan(numel(raw),1); % Create NaN column of appropriate size
            for ii = 1:numel(raw) 
                try
                    v(ii) = str2double(raw{ii});
                catch
                    v(ii) = NaN; % Assign NaN if individual conversion fails
                end
            end
            % warning('parseNumeric: Could not convert all cell elements to double.');
        end
    elseif ischar(raw) || isstring(raw) % Handle char array or string scalar
        parts = split(string(raw), ','); parts(parts=='') = [];
        try
            v = str2double(parts);
        catch
            v = nan(size(parts));
            % warning('parseNumeric: Could not convert all string parts to double.');
        end
    else
        v = []; % Unknown type
        % warning('parseNumeric: Unknown data type encountered.');
    end
    % Ensure it's a column vector
    if ~iscolumn(v) && ~isempty(v) % Add ~isempty(v) check
        v = v';
    end
end