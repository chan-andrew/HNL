function realtimestampingv2()

    % 1) prompt for JSON
    [fileName, filePath] = uigetfile('*.json', 'Select BrainSense JSON file');
    if isequal(fileName, 0)
        fprintf('No file selected. Exiting.\n');
        return;
    end
    fullFile = fullfile(filePath, fileName);
    try
        jsonText = fileread(fullFile);
        data = jsondecode(jsonText);
    catch ME
        error('Error reading or decoding JSON: %s', ME.message);
    end

    maxSeq = 255;

    % 2) process BrainSenseTimeDomain
    if isfield(data, 'BrainSenseTimeDomain')
        TD = data.BrainSenseTimeDomain;
        for k = 1:numel(TD)
            % anchor to UTC timestamp of first packet
            t0 = datetime(TD(k).FirstPacketDateTime, ...
                'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z''', 'TimeZone', 'UTC');

            % parse TicksInMses
            rawTicks = TD(k).TicksInMses;
            if isnumeric(rawTicks)
                ticks = rawTicks(:);
            elseif iscell(rawTicks)
                ticks = str2double(rawTicks(:));
            else
                parts = split(string(rawTicks), ',');
                parts(parts == "") = [];
                ticks = str2double(parts);
            end

            % parse GlobalPacketSizes
            rawSizes = TD(k).GlobalPacketSizes;
            if isnumeric(rawSizes)
                sizes = rawSizes(:);
            elseif iscell(rawSizes)
                sizes = str2double(rawSizes(:));
            else
                parts = split(string(rawSizes), ',');
                parts(parts == "") = [];
                sizes = str2double(parts);
            end

            % parse GlobalSequences
            rawSeqs = TD(k).GlobalSequences;
            if isnumeric(rawSeqs)
                seqs = rawSeqs(:);
            elseif iscell(rawSeqs)
                seqs = str2double(rawSeqs(:));
            else
                parts = split(string(rawSeqs), ',');
                parts(parts == "") = [];
                seqs = str2double(parts);
            end

            fs = TD(k).SampleRateInHz;
            dt = seconds(1/fs);

            % compute packet start times
            tStart = t0 + milliseconds(ticks - ticks(1));

            % detect sequence drops
            dseq = diff(seqs);
            dseq(dseq < 0) = dseq(dseq < 0) + maxSeq + 1;
            dropIdx = find(dseq ~= 1);
            if ~isempty(dropIdx)
                warning('TD(%d): Dropped packet(s) at indices %s', k, mat2str(dropIdx));
            end

            % detect tick misalignments
            dticks = diff(ticks);
            expectedGap = (sizes(1:end-1)/fs) * 1000;  % in ms
            misIdx = find(abs(dticks - expectedGap) > 2);
            if ~isempty(misIdx)
                warning('TD(%d): Misaligned packet(s) at indices %s', k, mat2str(misIdx));
            end

            % build timestamps per packet
            segs = cell(numel(ticks), 1);
            for i = 1:numel(ticks)
                if i > 1 && abs(dticks(i-1) - expectedGap(i-1)) <= 2
                    % continuous block
                    segs{i} = (segs{i-1}(end) + dt : dt : segs{i-1}(end) + dt * sizes(i)).';
                else
                    % re-anchor on packet start
                    segs{i} = (tStart(i) : dt : tStart(i) + dt * (sizes(i)-1)).';
                end
            end

            % concatenate and assign
            TD(k).TimeStamps = vertcat(segs{:});
            TD(k).TimeStamps.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
        end
        assignin('base', 'BrainSenseTimeDomain', TD);
    else
        warning('No field "BrainSenseTimeDomain" in JSON.');
    end

    % 3) process LFP streams if present
    lfpVars = {'BrainSenseLfp', 'LFPData', 'LFPmontageTimeDomain'};
    lfpName = '';
    for n = 1:numel(lfpVars)
        if isfield(data, lfpVars{n})
            lfpName = lfpVars{n};
            break;
        end
    end
    if isempty(lfpName)
        warning('No LFP struct found (BrainSenseLfp, LFPData, etc.) in JSON.');
    else
        LFP = data.(lfpName);
        for k = 1:numel(LFP)
            % anchor timestamp
            t0 = datetime(LFP(k).FirstPacketDateTime, ...
                'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z''', 'TimeZone', 'UTC');

            % parse LfpData.TicksInMs
            rawTicks = LFP(k).LfpData.TicksInMs;
            if isnumeric(rawTicks)
                ticks = rawTicks(:);
            elseif iscell(rawTicks)
                ticks = str2double(rawTicks(:));
            else
                parts = split(string(rawTicks), ',');
                parts(parts == "") = [];
                ticks = str2double(parts);
            end

            % parse GlobalPacketSizes or default to ones
            if isfield(LFP(k), 'GlobalPacketSizes')
                rawSizes = LFP(k).GlobalPacketSizes;
                if isnumeric(rawSizes)
                    sizes = rawSizes(:);
                elseif iscell(rawSizes)
                    sizes = str2double(rawSizes(:));
                else
                    parts = split(string(rawSizes), ',');
                    parts(parts == "") = [];
                    sizes = str2double(parts);
                end
            else
                sizes = ones(size(ticks));
            end

            % parse GlobalSequences if present and detect drops
            if isfield(LFP(k), 'GlobalSequences')
                rawSeqs = LFP(k).GlobalSequences;
                if isnumeric(rawSeqs)
                    seqs = rawSeqs(:);
                elseif iscell(rawSeqs)
                    seqs = str2double(rawSeqs(:));
                else
                    parts = split(string(rawSeqs), ',');
                    parts(parts == "") = [];
                    seqs = str2double(parts);
                end
                dseq = diff(seqs);
                dseq(dseq < 0) = dseq(dseq < 0) + maxSeq + 1;
                dropIdx = find(dseq ~= 1);
                if ~isempty(dropIdx)
                    warning('LFP(%d): Dropped packet(s) at indices %s', k, mat2str(dropIdx));
                end
            end

            fs = LFP(k).SampleRateInHz;
            dt = seconds(1/fs);
            tStart = t0 + milliseconds(ticks - ticks(1));

            % build timestamps
            segs = cell(numel(ticks), 1);
            for i = 1:numel(ticks)
                segs{i} = (tStart(i) : dt : tStart(i) + dt * (sizes(i)-1)).';
            end

            LFP(k).TimeStamps = vertcat(segs{:});
            LFP(k).TimeStamps.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
        end
        assignin('base', lfpName, LFP);
    end

    disp('All TimeStamps successfully assigned.');
end