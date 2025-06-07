function Chan_LoadJson()
    % Load a BrainSense JSON export and assign top-level fields to workspace

    % Ask user to select a .json file
    [file, path] = uigetfile('*.json', 'Select a BrainSense JSON file');
    if isequal(file, 0)
        fprintf('No file selected. Exiting.\n');
        return;
    end
    filePath = fullfile(path, file);

    % Read and decode JSON file
    try
        jsonText = fileread(filePath);
        data = jsondecode(jsonText);
    catch ME
        fprintf('Error reading or decoding file: %s\n', ME.message);
        return;
    end

    % Loop through each top-level field and assign it to base workspace
    fn = fieldnames(data);
    for i = 1:numel(fn)
        assignin('base', fn{i}, data.(fn{i}));
        fprintf('Assigned variable: %s\n', fn{i});
    end

    fprintf('\n Loaded %d variables into workspace from %s\n', numel(fn), file);
end
