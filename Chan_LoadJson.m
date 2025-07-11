function organizeJsonReports(srcFolder)

    % prompting user to select source and destination folder
    fprintf('Please select the folder containing the JSON files you want to rename.\n');
    srcFolder = uigetdir;
    fprintf('Please select the folder which you want the copied and renamed JSON files to go to.\n');
    destRoot = uigetdir;

    % using slightly modified LoadJson.m to load the json files
    [jsondat, fn] = LoadJson(srcFolder);
    if isempty(fn) % single file case
        fn = {srcFolder};
        jsondat = {jsondat};
        srcFolder = fileparts(srcFolder);
    end

    % processing
    for k = 1:numel(fn)
        src = fn{k};
        data = jsondat{k};

        % a) just pulling last name data
        lastName = data.PatientInformation.Initial.PatientLastName;

        % b) visit date: dd-mm-yyyy
        rawDate = data.SessionDate;
        dt = datetime(rawDate,'InputFormat','yyyy-MM-dd''T''HH:mm:ss''Z''', 'TimeZone','UTC');
        visitDate = datestr(dt, 'dd-mm-yyyy');
        %dt.Format = 'dd-MM-yyyy';
        %visitDate = string(dt);

        % c) diagnosis (also strips the prefix so there is no "DiagnosisTypeDef."
        rawDiag = data.PatientInformation.Initial.Diagnosis;  
        parts = split(rawDiag, '.');
        diagnosis = parts{end};

        % d) build base name & destination folder
        baseName = sprintf('%s_%s_%s', lastName, visitDate, diagnosis);
        destFolder = fullfile(destRoot, sprintf('%s_%s', lastName, diagnosis));
        if ~exist(destFolder, 'dir')
            mkdir(destFolder);
        end
        ext = '.json';
        destPath = fullfile(destFolder, [baseName ext]);

        % duplicate cleaning and collisions
        srcText = fileread(src);
        if exist(destPath, 'file')
            destText = fileread(destPath);
            if strcmp(srcText, destText)
                fprintf('Skipping duplicate (same content): %s\n', src);
                continue;
            else
                % collisions prevention if needed
                i = 1;
                while true
                    altName = sprintf('%s_%d%s', baseName, i, ext);
                    altPath = fullfile(destFolder, altName);
                    if ~exist(altPath, 'file')
                        destPath = altPath;
                        break;
                    end
                    if strcmp(srcText, fileread(altPath))
                        fprintf('Skipping duplicate in alt slot: %s\n', src);
                        break; 
                    end
                    i = i + 1;
                end
            end
        end

        % moving file to detination folder
        copyfile(src, destPath);
        fprintf('Copied:\n  %s\n→ %s\n', src, destPath);
    end

    fprintf('Completed.\n');

end
