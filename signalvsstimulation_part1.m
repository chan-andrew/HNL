function final = signalvsstimulation()

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

    % Initialize accumulators for plotting
    acc_orig_Lstim_times_utc = [];
    acc_orig_Lstim_values = [];
    acc_orig_Rstim_times_utc = [];
    acc_orig_Rstim_values = [];

    % 4) Align channels by FirstPacketDateTime (for warning only)
    tdFP  = datetime({TD.FirstPacketDateTime},  'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
    lfpFP = datetime({LFP.FirstPacketDateTime}, 'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
    idxMap = nan(size(lfpFP));
    for i = 1:numel(lfpFP)
        if ~isnat(lfpFP(i))
            [~, idxMap(i)] = min(abs(tdFP - lfpFP(i)));
        end
    end
    if numel(unique(idxMap(~isnan(idxMap)))) ~= numel(idxMap(~isnan(idxMap)))
        warning('One or more LFP channels align to the same TD channel based on FirstPacketDateTime proximity.');
    end

    % 5) Extract raw and stim per channel pair
    final = table();
    for j = 1:numel(LFP)
        % Determine channel pairing
        chanPair = split(LFP(j).Channel, ',');
        chanLeft = strtrim(chanPair{1});
        chanRight = '';
        if numel(chanPair) > 1, chanRight = strtrim(chanPair{2}); end

        fpStr = LFP(j).FirstPacketDateTime;
        idxL = find(strcmp({TD.Channel}, chanLeft) & strcmp({TD.FirstPacketDateTime}, fpStr));
        idxR = [];
        if ~isempty(chanRight)
            idxR = find(strcmp({TD.Channel}, chanRight) & strcmp({TD.FirstPacketDateTime}, fpStr));
        end
        if numel(idxL)~=1 || numel(idxR)~=1
            continue;
        end

        % Get time series and raw signals
        tvec = TD(idxL).TimeStamps;
        rawLeft = TD(idxL).TimeDomainData(:);
        rawRight = TD(idxR).TimeDomainData(:);

        % Initialize local stim arrays for this channel
        current_maL = [];
        current_maR = [];

        % Collect stimulation events for this LFP stream
        if isfield(LFP(j),'LfpData') && ~isempty(LFP(j).LfpData)
            pkt_ticks = arrayfun(@(pkt) pkt.TicksInMs, LFP(j).LfpData);
            t0_LFP = datetime(LFP(j).FirstPacketDateTime,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
            evt_times = t0_LFP + milliseconds(pkt_ticks - pkt_ticks(1));

            if isfield(LFP(j).LfpData(1),'Left') && isfield(LFP(j).LfpData(1).Left,'mA')
                current_maL = arrayfun(@(pkt) pkt.Left.mA, LFP(j).LfpData);
                acc_orig_Lstim_times_utc = [acc_orig_Lstim_times_utc; evt_times(:)];
                acc_orig_Lstim_values = [acc_orig_Lstim_values; current_maL(:)];
            end
            if isfield(LFP(j).LfpData(1),'Right') && isfield(LFP(j).LfpData(1).Right,'mA')
                current_maR = arrayfun(@(pkt) pkt.Right.mA, LFP(j).LfpData);
                acc_orig_Rstim_times_utc = [acc_orig_Rstim_times_utc; evt_times(:)];
                acc_orig_Rstim_values = [acc_orig_Rstim_values; current_maR(:)];
            end
        end

        % Expand stim vectors to match raw data length
        Nsig_L = numel(rawLeft);
        Nsig_R = numel(rawRight);
        stimL = nan(Nsig_L,1);
        stimR = nan(Nsig_R,1);
        if ~isempty(current_maL)
            repsL = ceil(Nsig_L/numel(current_maL));
            stim_fullL = repelem(current_maL, repsL);
            stimL = stim_fullL(1:Nsig_L);
        end
        if ~isempty(current_maR)
            repsR = ceil(Nsig_R/numel(current_maR));
            stim_fullR = repelem(current_maR, repsR);
            stimR = stim_fullR(1:Nsig_R);
        end

        % Append to final table
        T = table(tvec, rawLeft, stimL, rawRight, stimR, ...
                  'VariableNames',{'Time','LeftRaw','LeftStim_mA','RightRaw','RightStim_mA'});
        final = [final; T];
    end

    % 6) Assign to workspace and sort
    assignin('base','final',final);
    if ~isempty(final)
        disp('''final'' table created in workspace.');
        if any(diff(final.Time) < 0)
            final = sortrows(final,'Time');
            disp('Sorted final table by time.');
        end
    else
        disp('''final'' table is empty. Check LFP/TD matching.');
    end

    % 7) Plot coregistered TD & stim intensity
    figure('Name','Signal vs Stimulation');
    subplot(2,1,1);
      yyaxis left;
      plot(final.LeftRaw);
      ylabel('Left Raw Signal');
      yyaxis right;
      plot(final.LeftStim_mA);
      ylabel('Left Stim (mA)');
      title('Left Signal & Stim Intensity');
      xlabel('Time'); grid on;
    subplot(2,1,2);
      yyaxis left;
      plot(final.RightRaw);
      ylabel('Right Raw Signal');
      yyaxis right;
      plot(final.RightStim_mA);
      ylabel('Right Stim (mA)');
      title('Right Signal & Stim Intensity');
      xlabel('Time'); grid on;
end

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

        fs = S(k).SampleRateInHz;
        if fs == 0 % Avoid division by zero if SampleRateInHz is 0
            warning('SampleRateInHz is 0 for %s channel %d (%s). Timestamps may be incorrect.', type, k, S(k).Channel);
            dt = seconds(0); % Or handle as error / assign NaN
        else
            dt = seconds(1/fs);
        end
        
        t0 = datetime(S(k).FirstPacketDateTime, ...
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
