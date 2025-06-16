clear

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

% 6) Loop & match ±1.5 s
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
    
    % grab the correct side’s RateInHertz
    snap = BrainSenseLfp(idx).TherapySnapshot;
    if endsWith(tdChan,'_LEFT')
        rate = snap.Left.RateInHertz;
    else
        rate = snap.Right.RateInHertz;
    end
    
    BrainSenseTimeDomain.StimRateHz(i) = rate;
end

% 8) Build a “ThreePart” column of openable sub‐tables containing DateTime (placeholders for StimAmp & GammaPower)
nTD = height(BrainSenseTimeDomain);
% convert the original ISO strings to datetime once
DateTime = datetime(BrainSenseTimeDomain.FirstPacketDateTime, ...
    'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');

% preallocate cell array of subtables
threePartCells = cell(nTD,1);
for i = 1:nTD
    % each cell gets a 1×3 table:
    %   DateTime – actual timestamp
    %   StimAmp – NaN placeholder
    %   GammaPower – NaN placeholder
    threePartCells{i} = table( ...
        DateTime(i), ...
        nan, ...  % fill in later with stimulation amplitude (mA)
        nan, ...  % fill in later with gamma‐band power
        'VariableNames', {'DateTime','StimAmp','GammaPower'} );
end

% attach it as a new column to your main table
BrainSenseTimeDomain.ThreePart = threePartCells;

% 9) Compute average stimulation amplitude (mA) into ThreePart.StimAmp
for i = 1:nTD
    % re‐match the LFP packet for this time‐domain row
    tdChan = BrainSenseTimeDomain.Channel{i};
    tdTime = tdTimes(i);
    dt  = abs(lfpTimes - tdTime);
    idx = find(dt <= tol, 1);
    if isempty(idx)
        continue;
    end
    
    % make sure the side matches
    lfpChs = strsplit(BrainSenseLfp(idx).Channel, ',');
    if ~any(strcmp(lfpChs, tdChan))
        continue;
    end
    
    % pull out the mA readings
    L = BrainSenseLfp(idx).LfpData;  % M×1 struct array
    if endsWith(tdChan,'_LEFT')
        vals = arrayfun(@(s) s.Left.mA,  L);
    else
        vals = arrayfun(@(s) s.Right.mA, L);
    end
    
    % average and store
    avg_mA = mean(vals);
    tmpTable = BrainSenseTimeDomain.ThreePart{i};
    tmpTable.StimAmp = avg_mA;
    BrainSenseTimeDomain.ThreePart{i} = tmpTable;
end

% 10) Push updated table back to base workspace
assignin('base','BrainSenseTimeDomain',BrainSenseTimeDomain);

fprintf('Done.\n');