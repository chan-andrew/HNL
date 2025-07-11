function organizeJsonReports(srcFolder)

    % prompting user to select source and destination folder
    fprintf('Please select the folder containing the JSON files you want to rename.\n');
    srcFolder = uigetdir;
    fprintf('Please select the folder which you want the copied and renamed JSON files to go to.\n');
    destRoot = uigetdir;

    % Recursively get list of JSON files from all subdirectories
    jsonFiles = getAllJsonFiles(srcFolder);
    if isempty(jsonFiles)
        fprintf('No JSON files found in the selected folder or its subdirectories.\n');
        return;
    end
    
    fprintf('Found %d JSON files. Processing...\n', length(jsonFiles));
    
    % Process each file individually with error handling
    successCount = 0;
    errorCount = 0;
    
    for i = 1:length(jsonFiles)
        src = jsonFiles(i).fullPath;  % Use full path instead of constructing it
        fprintf('Processing file %d/%d: %s\n', i, length(jsonFiles), jsonFiles(i).name);
        
        try
            % Try to read and parse the JSON file
            txt = fileread(src);
            data = jsondecode(txt);
            
            % Validate required fields exist
            if ~isfield(data, 'PatientInformation') || ...
               ~isfield(data.PatientInformation, 'Initial') || ...
               ~isfield(data.PatientInformation.Initial, 'PatientLastName') || ...
               ~isfield(data.PatientInformation.Initial, 'Diagnosis') || ...
               ~isfield(data, 'SessionDate')
                fprintf('  WARNING: Missing required fields in %s, skipping...\n', jsonFiles(i).name);
                errorCount = errorCount + 1;
                continue;
            end
            
            % a) just pulling last name data
            lastName = data.PatientInformation.Initial.PatientLastName;

            % b) visit date: dd-mm-yyyy
            rawDate = data.SessionDate;
            dt = datetime(rawDate,'InputFormat','yyyy-MM-dd''T''HH:mm:ss''Z''', 'TimeZone','UTC');
            visitDate = datestr(dt, 'dd-mm-yyyy');

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
                    fprintf('  Skipping duplicate (same content): %s\n', src);
                    continue;
                else
                    % collisions prevention if needed
                    j = 1;
                    while true
                        altName = sprintf('%s_%d%s', baseName, j, ext);
                        altPath = fullfile(destFolder, altName);
                        if ~exist(altPath, 'file')
                            destPath = altPath;
                            break;
                        end
                        if strcmp(srcText, fileread(altPath))
                            fprintf('  Skipping duplicate in alt slot: %s\n', src);
                            break; 
                        end
                        j = j + 1;
                    end
                end
            end

            % moving file to destination folder
            copyfile(src, destPath);
            fprintf('  Copied:\n    %s\n  → %s\n', src, destPath);
            successCount = successCount + 1;
            
        catch ME
            fprintf('  ERROR processing %s: %s\n', jsonFiles(i).name, ME.message);
            errorCount = errorCount + 1;
            continue;
        end
    end

    fprintf('\nCompleted.\n');
    fprintf('Successfully processed: %d files\n', successCount);
    fprintf('Errors encountered: %d files\n', errorCount);

end

function jsonFiles = getAllJsonFiles(rootFolder)
    % Recursively find all JSON files in rootFolder and all subdirectories
    jsonFiles = [];
    
    % Get all items in current folder
    items = dir(rootFolder);
    
    for i = 1:length(items)
        item = items(i);
        
        % Skip . and .. directories
        if strcmp(item.name, '.') || strcmp(item.name, '..')
            continue;
        end
        
        fullPath = fullfile(rootFolder, item.name);
        
        if item.isdir
            % If it's a directory, recursively search it
            subJsonFiles = getAllJsonFiles(fullPath);
            jsonFiles = [jsonFiles, subJsonFiles];
        elseif endsWith(item.name, '.json', 'IgnoreCase', true)
            % If it's a JSON file, add it to our list
            jsonFileStruct = struct();
            jsonFileStruct.name = item.name;
            jsonFileStruct.fullPath = fullPath;
            jsonFileStruct.folder = rootFolder;
            jsonFiles = [jsonFiles, jsonFileStruct];
        end
    end

end
