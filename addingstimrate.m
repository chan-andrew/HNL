% 1) Prompt for JSON file and decode
[fileName, filePath] = uigetfile('*.json','Select a BrainSense JSON file');
if isequal(fileName,0)
    error('No file selected. Exiting.');
end
jsonText = fileread(fullfile(filePath,fileName));
data     = jsondecode(jsonText);

% 2) Extract the two fields
if ~isfield(data,'BrainSenseTimeDomain') || ~isfield(data,'BrainSenseLfp')
    error('JSON missing required fields.');
end

rawTD = data.BrainSenseTimeDomain;   % often a struct array
BrainSenseLfp = data.BrainSenseLfp;  % struct array

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
tdTimes  = datetime(BrainSenseTimeDomain.FirstPacketDateTime, ...
             'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
nLFP     = numel(BrainSenseLfp);
lfpTimes = datetime({BrainSenseLfp.FirstPacketDateTime}, ...
             'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');

% 6) Loop & match ±1.5 s
tol = seconds(1.5);
for i = 1:nTD
    tdChan = BrainSenseTimeDomain.Channel{i};
    tdTime = tdTimes(i);
    
    % find the first LFP row within tolerance
    dt  = abs(lfpTimes - tdTime);
    idx = find(dt <= tol, 1);
    if isempty(idx)
        continue;
    end
    
    % split LFP channels, pick the one matching tdChan
    lfpChs = strsplit(BrainSenseLfp(idx).Channel, ',');
    if ~any(strcmp(lfpChs, tdChan))
        continue;
    end
    
    % grab the correct side’s RateInHertz
    snap = BrainSenseLfp(idx).TherapySnapshot;
    if endsWith(tdChan,'_LEFT')
        rate = snap.Left.RateInHertz;
    else
        rate = snap.Right.RateInHertz;
    end
    
    BrainSenseTimeDomain.StimRateHz(i) = rate;
end

% 7) Push back to base
assignin('base','BrainSenseTimeDomain',BrainSenseTimeDomain);

fprintf('Done—BrainSenseTimeDomain now has a StimRateHz column in your workspace.\n');
