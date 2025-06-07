% [timestampsAll, dropFlagsAll] = droppedpacketdetectionworking(BrainSenseTimeDomain, true);
% tests for dopped packets

function [timestampsAll, dropFlagsAll] = droppedpacketdetectionworking(TDArray, doPlot)
    if nargin<2 || isempty(doPlot)
        doPlot = false;
    end
    nCh = numel(TDArray);
    timestampsAll = cell(nCh,1);
    dropFlagsAll = cell(nCh,1);

    % process each channel/field
    for ch = 1:nCh
        TD = TDArray(ch);
        [ts, drops] = assignSingle(TD);
        timestampsAll{ch} = ts;
        dropFlagsAll{ch} = drops;
    end

    % plot packet-drop flags per channel
    if doPlot
        figure;
        nCols = ceil(sqrt(nCh));
        nRows = ceil(nCh / nCols);
        for ch = 1:nCh
            subplot(nRows, nCols, ch);
            drops = dropFlagsAll{ch};
            stem(1:numel(drops), double(drops), 'filled');
            xlabel('Packet Index'); ylabel('Dropped (1=true)');
            title(sprintf('Channel %d', ch));
            ylim([-0.1,1.1]); grid on;
        end
        sgtitle('Packet Drop Flags per Channel');
    end
end

% sub-function: handle single struct
function [timestamps, dropPackets] = assignSingle(TD)
    % parse sizes
    rawSizes = TD.GlobalPacketSizes;
    if ischar(rawSizes) || isstring(rawSizes)
        sizes = str2double(strsplit(char(rawSizes), ','));
    else
        sizes = double(rawSizes(:)');
    end
    % parse ticks
    rawTicks = TD.TicksInMses;
    if ischar(rawTicks) || isstring(rawTicks)
        ticks = str2double(strsplit(char(rawTicks), ','));
    else
        ticks = double(rawTicks(:)');
    end
    fs = double(TD.SampleRateInHz);
    dt_ms = 1000 / fs;

    nPackets = numel(ticks);
    dropPackets = false(nPackets,1);
    time_ms = [];

    % build relative time (ms) backward
    for i = nPackets:-1:1
        N = sizes(i);
        tLast = ticks(i);
        if isempty(time_ms)
            t = (tLast - (N-1)*dt_ms) : dt_ms : tLast;
        else
            nextTick = ticks(i+1);
            if (nextTick - tLast) > ((N+1)*dt_ms)
                dropPackets(i) = true;
                t = (tLast - (N-1)*dt_ms) : dt_ms : tLast;
            else
                t_end = time_ms(1) - dt_ms;
                t_start = t_end - (N-1)*dt_ms;
                t = t_start : dt_ms : t_end;
            end
        end
        time_ms = [t, time_ms];
    end

    % Convert to absolute timestamps
    t0 = datetime(TD.FirstPacketDateTime, 'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''', 'TimeZone','UTC');
    rel_ms = time_ms - ticks(1);
    timestamps = t0 + milliseconds(rel_ms);
end
