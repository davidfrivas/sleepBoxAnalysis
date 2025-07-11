% This script:
% 1. Finds all CSV files in input folder (no subfolders)
% 2. Groups mice by user-defined genotypes (can be more than 2)
% 3. Retrieves epoch durations for each mouse from CSV files
% 4. Sorts sleep bouts into light phase (6:00-18:00) and dark phase (18:00-6:00)
% 5. Bins the epochs by their duration (2s, 4s, 8s, 16s, 32s, 64s, 128s, 256s, 512s)
% 6. Performs ZT (Zeitgeber Time) hour analysis, where ZT0 = 6:00, ZT12 = 18:00
% 7. SEPARATES DATA BY EXPERIMENTAL DAY (Day 1, Day 2, etc.) - days run from 6AM to 6AM
% 8. Computes averages and standard deviations by genotype (overall, light phase, dark phase, ZT hour, per day)
% 9. Plots results as dot graphs with error bars for all data categories
% 10. EXPORTS ALL DATA TO EXCEL WITH MULTIPLE SHEETS

%% Helper function to create valid MATLAB field names
function validName = makeValidFieldName(inputName)
    % Replace hyphens with underscores and ensure it starts with a letter
    validName = strrep(inputName, '-', '_');
    % If it starts with a number, prepend 'mouse_'
    if ~isempty(validName) && isstrprop(validName(1), 'digit')
        validName = ['mouse_' validName];
    end
    % Remove any other invalid characters (keep only letters, numbers, underscores)
    validName = regexprep(validName, '[^a-zA-Z0-9_]', '_');
end

%% Parameters and Setup
mainFolder = '/Users/davidrivas/Documents/research/sleep-box/Chd8Het_Old_July2023/091823'; % Use current directory, or specify your path

% Define bin edges for epoch durations (in seconds)
binEdges = [2, 4, 8, 16, 32, 64, 128, 256, 512, inf];
binLabels = {'2-4s', '4-8s', '8-16s', '16-32s', '32-64s', '64-128s', '128-256s', '256-512s', '>512s'};
numBins = length(binLabels);

% Create output folder for saving figures
outputFolder = fullfile(mainFolder, 'compiled_plots');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
    fprintf('Created output folder: %s\n', outputFolder);
else
    fprintf('Using existing output folder: %s\n', outputFolder);
end

% Keep track of all days found across all mice
allDays = [];

%% Find all CSV files and extract animal IDs
csvFiles = dir(fullfile(mainFolder, '*SB_2sec_*.csv'));

if isempty(csvFiles)
    error('No CSV files found in the specified folder. Files should follow the pattern: MM-DD-YYSB_2sec_mouseID.csv');
end

% Extract animal IDs from file names
mouseIDs = cell(length(csvFiles), 1);
for i = 1:length(csvFiles)
    fileName = csvFiles(i).name;
    % Extract mouse ID after 'SB_2sec_'
    parts = split(fileName, 'SB_2sec_');
    if length(parts) > 1
        % Remove .csv extension
        mouseID = strrep(parts{2}, '.csv', '');
        mouseIDs{i} = mouseID;
    else
        % If the file doesn't follow expected naming convention
        fprintf('Warning: Could not extract mouse ID from filename: %s\n', fileName);
        mouseIDs{i} = ['Unknown_' num2str(i)];
    end
end

% Display all discovered animal IDs
fprintf('Discovered animal IDs:\n');
for i = 1:length(mouseIDs)
    fprintf('  %d. %s (from file: %s)\n', i, mouseIDs{i}, csvFiles(i).name);
end
fprintf('\n');

%% Prompt user to define genotypes
fprintf('Please enter the names of genotypes you want to analyze, separated by commas:\n');
fprintf('(Examples: wild-type, mutant OR control, treatment1, treatment2)\n');
genotypeInput = input('Genotypes: ', 's');

% Process the input - split by commas and trim whitespace
genotypeNames = strtrim(split(genotypeInput, ','));
for i = 1:length(genotypeNames)
    genotypeNames{i} = strtrim(genotypeNames{i});
end

% Remove any empty entries
genotypeNames = genotypeNames(~cellfun('isempty', genotypeNames));

if isempty(genotypeNames)
    error('No genotypes specified. Please restart the script and provide at least one genotype.');
end

numGenotypes = length(genotypeNames);
fprintf('\nYou have defined %d genotypes: %s\n\n', numGenotypes, strjoin(genotypeNames, ', '));

% Create field names for genotypes (for MATLAB struct compatibility)
genotypeFieldNames = cell(numGenotypes, 1);
for i = 1:numGenotypes
    genotypeFieldNames{i} = makeValidFieldName(genotypeNames{i});
end

%% Initialize data structures dynamically based on genotypes
% Create structures to store binned epoch counts by genotype
epochCounts = struct();
totalCounts = struct();
animalIDs = struct();
dayData = struct();
dayStats = struct();

for i = 1:numGenotypes
    genField = genotypeFieldNames{i};
    epochCounts.(genField) = struct('sleep', zeros(0, numBins), 'light', zeros(0, numBins), 'dark', zeros(0, numBins));
    totalCounts.(genField) = struct('sleep', [], 'light', [], 'dark', [], 'mouseIDs', {{}});
    animalIDs.(genField) = {};
    dayData.(genField) = struct();
    dayStats.(genField) = struct();
end

%% Collect mouse assignments for each genotype
genotypeAssignments = struct();
allAssignedMice = {};

for i = 1:numGenotypes
    genotypeName = genotypeNames{i};
    genField = genotypeFieldNames{i};
    
    fprintf('Please enter the IDs of %s mice, separated by commas:\n', genotypeName);
    miceInput = input(sprintf('%s mice: ', genotypeName), 's');
    
    % Process the input
    miceList = strtrim(split(miceInput, ','));
    for j = 1:length(miceList)
        miceList{j} = strtrim(miceList{j});
    end
    
    % Remove empty entries
    miceList = miceList(~cellfun('isempty', miceList));
    
    genotypeAssignments.(genField) = miceList;
    allAssignedMice = [allAssignedMice; miceList];
    
    fprintf('  %s: %s\n', genotypeName, strjoin(miceList, ', '));
end

% Verify inputs and warn about any unlisted mice
unlistedMice = setdiff(mouseIDs, allAssignedMice);
if ~isempty(unlistedMice)
    fprintf('\nWarning: The following mice were not assigned to any genotype and will be skipped:\n');
    for i = 1:length(unlistedMice)
        fprintf('  %s\n', unlistedMice{i});
    end
    fprintf('\n');
end

% Check for overlapping assignments
for i = 1:numGenotypes-1
    for j = i+1:numGenotypes
        genField1 = genotypeFieldNames{i};
        genField2 = genotypeFieldNames{j};
        overlap = intersect(genotypeAssignments.(genField1), genotypeAssignments.(genField2));
        if ~isempty(overlap)
            fprintf('\nERROR: The following mice were assigned to multiple genotypes (%s and %s):\n', ...
                genotypeNames{i}, genotypeNames{j});
            for k = 1:length(overlap)
                fprintf('  %s\n', overlap{k});
            end
            error('Mice cannot be assigned to multiple genotypes. Please restart the script.');
        end
    end
end

fprintf('\nProcessing will begin with:\n');
for i = 1:numGenotypes
    genField = genotypeFieldNames{i};
    fprintf('  %s: %s\n', genotypeNames{i}, strjoin(genotypeAssignments.(genField), ', '));
end
fprintf('\n');

%% Process each CSV file
for i = 1:length(csvFiles)
    csvPath = fullfile(csvFiles(i).folder, csvFiles(i).name);
    currentMouseID = mouseIDs{i};
    
    % Determine genotype based on user input
    mouseGenotype = '';
    genotypeDisplayName = '';
    for j = 1:numGenotypes
        genField = genotypeFieldNames{j};
        if ismember(currentMouseID, genotypeAssignments.(genField))
            mouseGenotype = genField;
            genotypeDisplayName = genotypeNames{j};
            break;
        end
    end
    
    if isempty(mouseGenotype)
        fprintf('Skipping mouse %s (no genotype assigned)...\n', currentMouseID);
        continue;
    end
    
    fprintf('Processing mouse %s (genotype: %s)...\n', currentMouseID, genotypeDisplayName);
    
    try
        % Read the CSV file with all necessary columns
        opts = detectImportOptions(csvPath);
        
        % Make sure we have at least 3 columns (timestamp, linear time, sleep bout duration)
        if size(opts.VariableNames, 2) >= 3
            % Select timestamp (first column) and sleep bout duration (third column)
            opts.SelectedVariableNames = opts.VariableNames([1, 3]);
            
            % Set variable types appropriately
            opts = setvartype(opts, opts.SelectedVariableNames{1}, 'string'); % timestamp as string
            opts = setvartype(opts, opts.SelectedVariableNames{2}, 'double'); % duration as double
            
            % Read the data
            sleepData = readtable(csvPath, opts);
            
            if ~isempty(sleepData)
                % Extract sleep bout durations
                sleepDurations = sleepData.(opts.SelectedVariableNames{2});
                % Remove NaN values
                validIdx = ~isnan(sleepDurations);
                sleepDurations = sleepDurations(validIdx);
                
                % Extract timestamps for valid durations
                timestamps = sleepData.(opts.SelectedVariableNames{1})(validIdx);
                
                % Parse timestamps and determine light/dark phase
                lightDurations = [];
                darkDurations = [];
                
                % Initialize day-specific storage for this mouse
                mouseDayData = struct();
                
                for t = 1:length(timestamps)
                    try
                        % Parse timestamp format (e.g., "9/19/23 6:01" or "2023-09-19 06:03:16")
                        timestamp = timestamps(t);
                        
                        % Parse the full datetime to determine experimental day
                        try
                            % Your timestamps are in format "2023-09-19 06:07:00"
                            timestampStr = char(timestamp);
                            
                            % Parse date and time manually for better compatibility
                            dateParts = split(timestampStr, ' ');
                            datePart = dateParts{1}; % "2023-09-19"
                            timePart = dateParts{2}; % "06:07:00"
                            
                            % Extract hour from time part
                            timeParts = split(timePart, ':');
                            hourNum = str2double(timeParts{1});
                            
                            % Convert to datenum for easier manipulation
                            fullDatenum = datenum(timestampStr, 'yyyy-mm-dd HH:MM:SS');
                            
                            % Calculate experimental day using datenum arithmetic
                            % Day boundaries are at 6:00 AM
                            if hourNum < 6
                                % Before 6 AM - belongs to previous day's experiment
                                expDayNum = floor(fullDatenum) - 1;
                            else
                                % 6 AM or later - belongs to current day's experiment
                                expDayNum = floor(fullDatenum);
                            end
                            
                            % Debug output for first few timestamps
                            if t <= 3
                                fprintf('    Debug: timestamp=%s, hour=%d, expDayNum=%.0f\n', ...
                                    timestampStr, hourNum, expDayNum);
                            end
                            
                        catch parseErr
                            fprintf('    Warning: Could not parse timestamp %s: %s\n', timestamp, parseErr.message);
                            continue;
                        end
                        
                        % Convert to day number (relative to first day found)
                        try
                            if isempty(allDays)
                                allDays = expDayNum;
                                dayNum = 1;
                                if t <= 3
                                    fprintf('    Debug: First day, dayNum=%d\n', dayNum);
                                end
                            else
                                % Find matching day using exact comparison (since we're using floor)
                                dayIdx = find(allDays == expDayNum);
                                
                                if isempty(dayIdx)
                                    allDays = [allDays; expDayNum];
                                    allDays = sort(allDays);
                                    dayIdx = find(allDays == expDayNum);
                                    if t <= 3
                                        fprintf('    Debug: New day found, dayIdx=%d\n', dayIdx(1));
                                    end
                                end
                                dayNum = dayIdx(1); % Take first match
                            end
                        catch dayErr
                            fprintf('    Warning: Error finding day number: %s\n', dayErr.message);
                            continue;
                        end
                        
                        % Initialize this day's data if needed
                        dayKey = sprintf('day%d', dayNum);
                        if ~isfield(mouseDayData, dayKey)
                            mouseDayData.(dayKey) = struct('sleep', [], 'light', [], 'dark', []);
                        end
                        
                        % Add this bout to the appropriate day
                        mouseDayData.(dayKey).sleep = [mouseDayData.(dayKey).sleep; sleepDurations(t)];
                        
                        % Extract hour for light/dark phase determination
                        timeParts = split(timestamp, ' ');
                        if length(timeParts) >= 2
                            timeStr = timeParts{2};
                            hourParts = split(timeStr, ':');
                            hour_orig = str2double(hourParts{1});
                            
                            % Determine phase based on hour (6:00-18:00 is light phase)
                            if hour_orig >= 6 && hour_orig < 18
                                % Light phase
                                lightDurations = [lightDurations; sleepDurations(t)];
                                mouseDayData.(dayKey).light = [mouseDayData.(dayKey).light; sleepDurations(t)];
                            else
                                % Dark phase
                                darkDurations = [darkDurations; sleepDurations(t)];
                                mouseDayData.(dayKey).dark = [mouseDayData.(dayKey).dark; sleepDurations(t)];
                            end
                        else
                            fprintf('    Warning: Unable to parse timestamp: %s\n', timestamp);
                        end
                    catch e
                        fprintf('    Warning: Error processing timestamp %s: %s\n', timestamp, e.message);
                    end
                end
                
                % Store the per-day data for this mouse
                % Convert mouse ID to valid MATLAB field name
                mouseFieldName = makeValidFieldName(currentMouseID);
                dayData.(mouseGenotype).(mouseFieldName) = mouseDayData;
                
                fprintf('  Found %d sleep bouts (%d light phase, %d dark phase)\n', ...
                    length(sleepDurations), length(lightDurations), length(darkDurations));
                
                % Bin epochs by duration for each category
                sleepBinCounts = histcounts(sleepDurations, binEdges);
                lightBinCounts = histcounts(lightDurations, binEdges);
                darkBinCounts = histcounts(darkDurations, binEdges);
                
                % Store the binned data for this mouse
                epochCounts.(mouseGenotype).sleep = [epochCounts.(mouseGenotype).sleep; sleepBinCounts];
                epochCounts.(mouseGenotype).light = [epochCounts.(mouseGenotype).light; lightBinCounts];
                epochCounts.(mouseGenotype).dark = [epochCounts.(mouseGenotype).dark; darkBinCounts];
                
                % Store total counts for this mouse
                totalCounts.(mouseGenotype).mouseIDs{end+1} = currentMouseID;
                totalCounts.(mouseGenotype).sleep(end+1) = length(sleepDurations);
                totalCounts.(mouseGenotype).light(end+1) = length(lightDurations);
                totalCounts.(mouseGenotype).dark(end+1) = length(darkDurations);
                
                % Store animal ID
                animalIDs.(mouseGenotype){end+1} = currentMouseID;
                
                % Display summary for this mouse
                fprintf('  Summary for mouse %s:\n', currentMouseID);
                fprintf('    All sleep bouts: %d (binned: %s)\n', length(sleepDurations), mat2str(sleepBinCounts));
                fprintf('    Light phase: %d (binned: %s)\n', length(lightDurations), mat2str(lightBinCounts));
                fprintf('    Dark phase: %d (binned: %s)\n', length(darkDurations), mat2str(darkBinCounts));
                
                % Display per-day summary
                dayFields = fieldnames(mouseDayData);
                fprintf('    Data found for %d days:\n', length(dayFields));
                for d = 1:length(dayFields)
                    dayKey = dayFields{d};
                    dayNumber = str2double(dayKey(4:end));
                    fprintf('      Day %d: %d sleep bouts (%d light, %d dark)\n', ...
                        dayNumber, ...
                        length(mouseDayData.(dayKey).sleep), ...
                        length(mouseDayData.(dayKey).light), ...
                        length(mouseDayData.(dayKey).dark));
                end
                
            else
                fprintf('  Warning: No sleep bout data found in file for mouse %s\n', currentMouseID);
            end
        else
            fprintf('  Warning: CSV file for mouse %s does not have the expected column structure\n', currentMouseID);
        end
    catch e
        fprintf('  Error processing file for mouse %s: %s\n', currentMouseID, e.message);
    end
end

%% Display overall day summary
fprintf('\n=== EXPERIMENTAL DAY SUMMARY ===\n');
fprintf('Total experimental days found: %d\n', length(allDays));
for d = 1:length(allDays)
    fprintf('  Day %d: %s\n', d, datestr(allDays(d), 'yyyy-mm-dd'));
end
fprintf('\n');

%% Calculate per-day statistics
fprintf('Calculating per-day statistics...\n');

numDays = length(allDays);

% For each genotype, calculate stats for each day
for i = 1:numGenotypes
    genField = genotypeFieldNames{i};
    
    % Initialize per-day storage
    for d = 1:numDays
        dayKey = sprintf('day%d', d);
        dayStats.(genField).(dayKey) = struct();
        dayStats.(genField).(dayKey).sleep = struct('counts', [], 'binned', []);
        dayStats.(genField).(dayKey).light = struct('counts', [], 'binned', []);
        dayStats.(genField).(dayKey).dark = struct('counts', [], 'binned', []);
    end
    
    % Process each mouse of this genotype
    if isfield(dayData, genField)
        mouseFields = fieldnames(dayData.(genField));
        for m = 1:length(mouseFields)
            mouseFieldName = mouseFields{m};
            mouseData = dayData.(genField).(mouseFieldName);
            
            % Process each day for this mouse
            mouseDayFields = fieldnames(mouseData);
            for d = 1:length(mouseDayFields)
                dayKey = mouseDayFields{d};
                
                % Get data for this day
                sleepBouts = mouseData.(dayKey).sleep;
                lightBouts = mouseData.(dayKey).light;
                darkBouts = mouseData.(dayKey).dark;
                
                % Store counts
                dayStats.(genField).(dayKey).sleep.counts = [dayStats.(genField).(dayKey).sleep.counts; length(sleepBouts)];
                dayStats.(genField).(dayKey).light.counts = [dayStats.(genField).(dayKey).light.counts; length(lightBouts)];
                dayStats.(genField).(dayKey).dark.counts = [dayStats.(genField).(dayKey).dark.counts; length(darkBouts)];
                
                % Bin the data
                sleepBinned = histcounts(sleepBouts, binEdges);
                lightBinned = histcounts(lightBouts, binEdges);
                darkBinned = histcounts(darkBouts, binEdges);
                
                % Store binned data
                dayStats.(genField).(dayKey).sleep.binned = [dayStats.(genField).(dayKey).sleep.binned; sleepBinned];
                dayStats.(genField).(dayKey).light.binned = [dayStats.(genField).(dayKey).light.binned; lightBinned];
                dayStats.(genField).(dayKey).dark.binned = [dayStats.(genField).(dayKey).dark.binned; darkBinned];
            end
        end
    end
    
    % Calculate means and SEMs for each day
    for d = 1:numDays
        dayKey = sprintf('day%d', d);
        
        % Sleep bouts
        if ~isempty(dayStats.(genField).(dayKey).sleep.counts)
            dayStats.(genField).(dayKey).sleep.mean_count = mean(dayStats.(genField).(dayKey).sleep.counts);
            dayStats.(genField).(dayKey).sleep.sem_count = std(dayStats.(genField).(dayKey).sleep.counts) / sqrt(length(dayStats.(genField).(dayKey).sleep.counts));
            dayStats.(genField).(dayKey).sleep.mean_binned = mean(dayStats.(genField).(dayKey).sleep.binned, 1);
            dayStats.(genField).(dayKey).sleep.sem_binned = std(dayStats.(genField).(dayKey).sleep.binned, 0, 1) / sqrt(size(dayStats.(genField).(dayKey).sleep.binned, 1));
        end
        
        % Light phase
        if ~isempty(dayStats.(genField).(dayKey).light.counts)
            dayStats.(genField).(dayKey).light.mean_count = mean(dayStats.(genField).(dayKey).light.counts);
            dayStats.(genField).(dayKey).light.sem_count = std(dayStats.(genField).(dayKey).light.counts) / sqrt(length(dayStats.(genField).(dayKey).light.counts));
            dayStats.(genField).(dayKey).light.mean_binned = mean(dayStats.(genField).(dayKey).light.binned, 1);
            dayStats.(genField).(dayKey).light.sem_binned = std(dayStats.(genField).(dayKey).light.binned, 0, 1) / sqrt(size(dayStats.(genField).(dayKey).light.binned, 1));
        end
        
        % Dark phase
        if ~isempty(dayStats.(genField).(dayKey).dark.counts)
            dayStats.(genField).(dayKey).dark.mean_count = mean(dayStats.(genField).(dayKey).dark.counts);
            dayStats.(genField).(dayKey).dark.sem_count = std(dayStats.(genField).(dayKey).dark.counts) / sqrt(length(dayStats.(genField).(dayKey).dark.counts));
            dayStats.(genField).(dayKey).dark.mean_binned = mean(dayStats.(genField).(dayKey).dark.binned, 1);
            dayStats.(genField).(dayKey).dark.sem_binned = std(dayStats.(genField).(dayKey).dark.binned, 0, 1) / sqrt(size(dayStats.(genField).(dayKey).dark.binned, 1));
        end
    end
end

%% Calculate statistics by genotype
% Initialize structures for mean and std values
meanCounts = struct();
stdCounts = struct();
semCounts = struct(); % Standard error of the mean

for i = 1:numGenotypes
    genField = genotypeFieldNames{i};
    
    if isfield(epochCounts, genField)
        % Calculate mean and std for all sleep bouts
        if ~isempty(epochCounts.(genField).sleep)
            meanCounts.(genField).sleep = mean(epochCounts.(genField).sleep, 1);
            stdCounts.(genField).sleep = std(epochCounts.(genField).sleep, 0, 1);
            semCounts.(genField).sleep = stdCounts.(genField).sleep / sqrt(size(epochCounts.(genField).sleep, 1));
        else
            fprintf('No sleep bout data found for genotype: %s\n', genotypeNames{i});
            meanCounts.(genField).sleep = zeros(1, numBins);
            stdCounts.(genField).sleep = zeros(1, numBins);
            semCounts.(genField).sleep = zeros(1, numBins);
        end
        
        % Calculate mean and std for light phase bouts
        if ~isempty(epochCounts.(genField).light)
            meanCounts.(genField).light = mean(epochCounts.(genField).light, 1);
            stdCounts.(genField).light = std(epochCounts.(genField).light, 0, 1);
            semCounts.(genField).light = stdCounts.(genField).light / sqrt(size(epochCounts.(genField).light, 1));
        else
            fprintf('No light phase data found for genotype: %s\n', genotypeNames{i});
            meanCounts.(genField).light = zeros(1, numBins);
            stdCounts.(genField).light = zeros(1, numBins);
            semCounts.(genField).light = zeros(1, numBins);
        end
        
        % Calculate mean and std for dark phase bouts
        if ~isempty(epochCounts.(genField).dark)
            meanCounts.(genField).dark = mean(epochCounts.(genField).dark, 1);
            stdCounts.(genField).dark = std(epochCounts.(genField).dark, 0, 1);
            semCounts.(genField).dark = stdCounts.(genField).dark / sqrt(size(epochCounts.(genField).dark, 1));
        else
            fprintf('No dark phase data found for genotype: %s\n', genotypeNames{i});
            meanCounts.(genField).dark = zeros(1, numBins);
            stdCounts.(genField).dark = zeros(1, numBins);
            semCounts.(genField).dark = zeros(1, numBins);
        end
    end
end

%% Export statistics to multiple sheets in an Excel file
fprintf('Creating Excel export with multiple sheets...\n');

% Create a consistent bin labels array that includes "Total"
allBinLabels = ['Total'; binLabels']; 

% 1. Create combined table for all sleep data (total + light/dark + genotype stats)
fprintf('Preparing combined total sleep data sheet...\n');
totalBinTable = array2table(zeros(length(allBinLabels), 0));
totalBinTable.Bin = allBinLabels;

% Add individual mouse data first
for i = 1:numGenotypes
    genField = genotypeFieldNames{i};
    
    if isfield(animalIDs, genField) && ~isempty(animalIDs.(genField))
        for m = 1:length(animalIDs.(genField))
            mouseID = animalIDs.(genField){m};
            fprintf('  Adding data for mouse: %s\n', mouseID);
            
            % Create data columns for all sleep phases
            if isfield(totalCounts, genField) && length(totalCounts.(genField).sleep) >= m
                sleepData = [totalCounts.(genField).sleep(m); epochCounts.(genField).sleep(m,:)'];
                lightData = [totalCounts.(genField).light(m); epochCounts.(genField).light(m,:)'];
                darkData = [totalCounts.(genField).dark(m); epochCounts.(genField).dark(m,:)'];
                
                % Add all columns for this mouse
                totalBinTable.([mouseID '_Sleep']) = sleepData;
                totalBinTable.([mouseID '_Light']) = lightData;
                totalBinTable.([mouseID '_Dark']) = darkData;
            end
        end
    end
end

% Add genotype statistics for all phases
for i = 1:numGenotypes
    genField = genotypeFieldNames{i};
    display_name = genotypeNames{i}; % Use original genotype name
    
    if isfield(totalCounts, genField)
        % Sleep statistics
        totalBinTable.([display_name '_Sleep_Mean']) = [mean(totalCounts.(genField).sleep); meanCounts.(genField).sleep'];
        totalBinTable.([display_name '_Sleep_SD']) = [std(totalCounts.(genField).sleep); stdCounts.(genField).sleep'];
        totalSEM_sleep = std(totalCounts.(genField).sleep)/sqrt(length(totalCounts.(genField).sleep));
        totalBinTable.([display_name '_Sleep_SEM']) = [totalSEM_sleep; semCounts.(genField).sleep'];
        
        % Light phase statistics
        totalBinTable.([display_name '_Light_Mean']) = [mean(totalCounts.(genField).light); meanCounts.(genField).light'];
        totalBinTable.([display_name '_Light_SD']) = [std(totalCounts.(genField).light); stdCounts.(genField).light'];
        totalSEM_light = std(totalCounts.(genField).light)/sqrt(length(totalCounts.(genField).light));
        totalBinTable.([display_name '_Light_SEM']) = [totalSEM_light; semCounts.(genField).light'];
        
        % Dark phase statistics
        totalBinTable.([display_name '_Dark_Mean']) = [mean(totalCounts.(genField).dark); meanCounts.(genField).dark'];
        totalBinTable.([display_name '_Dark_SD']) = [std(totalCounts.(genField).dark); stdCounts.(genField).dark'];
        totalSEM_dark = std(totalCounts.(genField).dark)/sqrt(length(totalCounts.(genField).dark));
        totalBinTable.([display_name '_Dark_SEM']) = [totalSEM_dark; semCounts.(genField).dark'];
    end
end

% Create per-day tables
fprintf('Preparing per-day data sheets...\n');

for d = 1:numDays
    dayKey = sprintf('day%d', d);
    
    % Create table for this day
    dayTable = array2table(zeros(length(allBinLabels), 0));
    dayTable.Bin = allBinLabels;
    
    % Add individual mouse data first
    fprintf('  Adding individual mouse data for Day %d...\n', d);
    for i = 1:numGenotypes
        genField = genotypeFieldNames{i};
        
        if isfield(animalIDs, genField) && ~isempty(animalIDs.(genField))
            for m = 1:length(animalIDs.(genField))
                mouseID = animalIDs.(genField){m};
                mouseFieldName = makeValidFieldName(mouseID);
                
                % Check if this mouse has data for this day
                if isfield(dayData, genField) && isfield(dayData.(genField), mouseFieldName)
                    mouseDayData = dayData.(genField).(mouseFieldName);
                    
                    if isfield(mouseDayData, dayKey)
                        % Get data for this specific day
                        sleepBouts = mouseDayData.(dayKey).sleep;
                        lightBouts = mouseDayData.(dayKey).light;
                        darkBouts = mouseDayData.(dayKey).dark;
                        
                        % Bin the data
                        sleepBinned = histcounts(sleepBouts, binEdges);
                        lightBinned = histcounts(lightBouts, binEdges);
                        darkBinned = histcounts(darkBouts, binEdges);
                        
                        % Create columns with total count + binned data
                        sleepData = [length(sleepBouts); sleepBinned'];
                        lightData = [length(lightBouts); lightBinned'];
                        darkData = [length(darkBouts); darkBinned'];
                        
                        % Add columns for this mouse
                        dayTable.([mouseID '_Sleep']) = sleepData;
                        dayTable.([mouseID '_Light']) = lightData;
                        dayTable.([mouseID '_Dark']) = darkData;
                        
                        fprintf('    Added data for mouse %s\n', mouseID);
                    else
                        % Mouse has no data for this day - add zeros
                        zeroData = zeros(length(allBinLabels), 1);
                        dayTable.([mouseID '_Sleep']) = zeroData;
                        dayTable.([mouseID '_Light']) = zeroData;
                        dayTable.([mouseID '_Dark']) = zeroData;
                        
                        fprintf('    Added zero data for mouse %s (no data for this day)\n', mouseID);
                    end
                else
                    % Mouse not found in dayData - add zeros
                    zeroData = zeros(length(allBinLabels), 1);
                    dayTable.([mouseID '_Sleep']) = zeroData;
                    dayTable.([mouseID '_Light']) = zeroData;
                    dayTable.([mouseID '_Dark']) = zeroData;
                    
                    fprintf('    Added zero data for mouse %s (not found in dayData)\n', mouseID);
                end
            end
        end
    end
    
    % Add genotype statistics (means and SEMs)
    fprintf('  Adding genotype statistics for Day %d...\n', d);
    for i = 1:numGenotypes
        genField = genotypeFieldNames{i};
        display_name = genotypeNames{i}; % Use original genotype name
        
        if isfield(dayStats, genField) && isfield(dayStats.(genField), dayKey)
            % Add mean and SEM for sleep, light, and dark phases
            phases = {'sleep', 'light', 'dark'};
            
            for p = 1:length(phases)
                phase = phases{p};
                
                if isfield(dayStats.(genField).(dayKey), phase) && ...
                   isfield(dayStats.(genField).(dayKey).(phase), 'mean_count')
                    
                    % Total count + binned means
                    meanData = [dayStats.(genField).(dayKey).(phase).mean_count; ...
                               dayStats.(genField).(dayKey).(phase).mean_binned'];
                    
                    % Total SEM + binned SEMs  
                    semData = [dayStats.(genField).(dayKey).(phase).sem_count; ...
                              dayStats.(genField).(dayKey).(phase).sem_binned'];
                    
                    dayTable.([display_name '_' phase '_Mean']) = meanData;
                    dayTable.([display_name '_' phase '_SEM']) = semData;
                end
            end
        end
    end
    
    % Write day table to Excel with new sheet name
    sheetName = sprintf('Day_%d_Bin', d);
    xlsFilePath = fullfile(mainFolder, 'sleep_bout_analysis.xls');
    writetable(dayTable, xlsFilePath, 'Sheet', sheetName);
    
    fprintf('  Added Day %d data to Excel sheet: %s\n', d, sheetName);
end

% 2. Create a summary table mapping mice to genotypes
mouseGenotypeSummary = table();
allMouseIDs_summary = {};
allGenotypes_summary = {};

for i = 1:numGenotypes
    genField = genotypeFieldNames{i};
    genotypeName = genotypeNames{i};
    
    if isfield(animalIDs, genField)
        allMouseIDs_summary = [allMouseIDs_summary; animalIDs.(genField)(:)];
        allGenotypes_summary = [allGenotypes_summary; repmat({genotypeName}, length(animalIDs.(genField)), 1)];
    end
end

mouseGenotypeSummary.MouseID = allMouseIDs_summary;
mouseGenotypeSummary.Genotype = allGenotypes_summary;

% Write all tables to different sheets in the Excel file
xlsFilePath = fullfile(mainFolder, 'sleep_bout_analysis.xls');
fprintf('Writing combined sheets to Excel file: %s\n', xlsFilePath);

% Write combined sheet - now includes all data (total + light/dark + genotype stats)
writetable(totalBinTable, xlsFilePath, 'Sheet', 'Total_Bin');
writetable(mouseGenotypeSummary, xlsFilePath, 'Sheet', 'Mouse_Genotype_Summary');

fprintf('Excel file created successfully.\n');

%% Create figures with dynamic colors and legends
% Turn off figure visibility to prevent windows from opening
set(0, 'DefaultFigureVisible', 'off');

% Generate distinct colors for each genotype
if numGenotypes <= 8
    % Use ColorBrewer-inspired colors for better visualization
    baseColors = [
        0.89, 0.10, 0.11;  % Red
        0.22, 0.49, 0.72;  % Blue  
        0.30, 0.69, 0.29;  % Green
        0.60, 0.31, 0.64;  % Purple
        1.00, 0.50, 0.00;  % Orange
        0.65, 0.34, 0.16;  % Brown
        0.97, 0.51, 0.75;  % Pink
        0.40, 0.40, 0.40;  % Gray
    ];
    genotypeColors = baseColors(1:numGenotypes, :);
else
    % Use lines colormap for many genotypes
    genotypeColors = colormap(lines(numGenotypes));
end

% Define the categories to plot
categories = {'sleep', 'light', 'dark'};
categoryLabels = {'All Sleep', 'Light Phase', 'Dark Phase'};

% Create figures for each category
for c = 1:length(categories)
    category = categories{c};
    categoryLabel = categoryLabels{c};
    
    %% Dot Plot (Individual mice with mean)
    figure('Name', [categoryLabel ' Bout Durations - Dot Plot'], 'Position', [100, 100, 800, 500]);
    hold on;
    
    % Initialize handles for legend
    plotHandles = [];
    legendTexts = {};
    
    % Plot data for each genotype
    for i = 1:numGenotypes
        genField = genotypeFieldNames{i};
        genotypeName = genotypeNames{i};
        color = genotypeColors(i, :);
        
        % Plot individual data points
        if isfield(epochCounts, genField) && ~isempty(epochCounts.(genField).(category))
            for j = 1:size(epochCounts.(genField).(category), 1)
                scatter((1:numBins) + (i-1)*0.1, epochCounts.(genField).(category)(j,:), 50, color, 'o', 'filled', 'MarkerFaceAlpha', 0.3);
            end
            
            % Plot mean with error bars and capture handle
            h = errorbar((1:numBins) + (i-1)*0.1, meanCounts.(genField).(category), semCounts.(genField).(category), ...
                'Color', color, 'LineWidth', 2, 'Marker', 'o', 'MarkerFaceColor', color, 'MarkerSize', 10);
            
            % Add to legend arrays
            plotHandles = [plotHandles, h];
            legendTexts{end+1} = sprintf('%s (Mean ± SEM)', genotypeName);
        end
    end
    
    % Add labels and legend
    title([categoryLabel ' Bout Durations']);
    xlabel('Bout Duration');
    ylabel('Number of Bouts');
    set(gca, 'XTick', 1:numBins, 'XTickLabel', binLabels);
    xtickangle(45);
    
    % Create legend with handles
    if ~isempty(plotHandles)
        legend(plotHandles, legendTexts, 'Location', 'best');
    end
    
    grid on;
    hold off;
    
    % Save the figure
    saveas(gcf, fullfile(outputFolder, [strrep(categoryLabel, ' ', '_') '_Bouts_Dot.fig']));
    saveas(gcf, fullfile(outputFolder, [strrep(categoryLabel, ' ', '_') '_Bouts_Dot.png']));
end

fprintf('All figures saved to %s\n', outputFolder);

%% Create per-day plots with dynamic genotypes
fprintf('Creating per-day plots...\n');

% Plot sleep bout counts across days
figure('Name', 'Sleep Bout Counts Across Days', 'Position', [100, 100, 1000, 600]);
hold on;

dayNumbers = 1:numDays;

% Initialize handles for legend
plotHandles = [];
legendTexts = {};

% Plot data for each genotype
for i = 1:numGenotypes
    genField = genotypeFieldNames{i};
    genotypeName = genotypeNames{i};
    color = genotypeColors(i, :);
    
    % Initialize arrays for this genotype
    genotypeMeans = zeros(1, numDays);
    genotypeSEMs = zeros(1, numDays);
    
    % Plot individual data points and collect means/SEMs
    for d = 1:numDays
        dayKey = sprintf('day%d', d);
        
        if isfield(dayStats.(genField), dayKey) && ~isempty(dayStats.(genField).(dayKey).sleep.counts)
            individualCounts = dayStats.(genField).(dayKey).sleep.counts;
            % Add small random jitter to x-position for visibility
            xPos = d + (rand(size(individualCounts)) - 0.5) * 0.1 + (i-1)*0.15;
            scatter(xPos, individualCounts, 50, color, 'o', 'filled', 'MarkerFaceAlpha', 0.3);
            
            % Store mean and SEM
            genotypeMeans(d) = dayStats.(genField).(dayKey).sleep.mean_count;
            genotypeSEMs(d) = dayStats.(genField).(dayKey).sleep.sem_count;
        end
    end
    
    % Plot means with error bars on top
    h = errorbar(dayNumbers + (i-1)*0.15, genotypeMeans, genotypeSEMs, 'Color', color, 'LineStyle', 'none', ...
        'Marker', 'o', 'MarkerFaceColor', color, 'MarkerSize', 10, 'LineWidth', 2, 'CapSize', 6);
    
    % Add to legend arrays
    plotHandles = [plotHandles, h];
    legendTexts{end+1} = sprintf('%s (Mean ± SEM)', genotypeName);
end

title('Sleep Bout Counts Across Experimental Days');
xlabel('Experimental Day');
ylabel('Number of Sleep Bouts');
legend(plotHandles, legendTexts, 'Location', 'best');
grid on;
hold off;

% Save figure
saveas(gcf, fullfile(outputFolder, 'Sleep_Bout_Counts_Across_Days.fig'));
saveas(gcf, fullfile(outputFolder, 'Sleep_Bout_Counts_Across_Days.png'));

% Create separate plots for light and dark phases across days
phases = {'light', 'dark'};
phaseLabels = {'Light Phase', 'Dark Phase'};

for p = 1:length(phases)
    phase = phases{p};
    phaseLabel = phaseLabels{p};
    
    figure('Name', [phaseLabel ' Sleep Bout Counts Across Days'], 'Position', [100, 100, 1000, 600]);
    hold on;
    
    % Initialize handles for legend
    plotHandles_phase = [];
    legendTexts_phase = {};
    
    % Plot data for each genotype
    for i = 1:numGenotypes
        genField = genotypeFieldNames{i};
        genotypeName = genotypeNames{i};
        color = genotypeColors(i, :);
        
        % Initialize arrays for this genotype and phase
        genotypeMeans_phase = zeros(1, numDays);
        genotypeSEMs_phase = zeros(1, numDays);
        
        % Plot individual data points and collect means/SEMs
        for d = 1:numDays
            dayKey = sprintf('day%d', d);
            
            if isfield(dayStats.(genField), dayKey) && isfield(dayStats.(genField).(dayKey).(phase), 'counts') && ...
               ~isempty(dayStats.(genField).(dayKey).(phase).counts)
                individualCounts = dayStats.(genField).(dayKey).(phase).counts;
                % Add small random jitter to x-position for visibility
                xPos = d + (rand(size(individualCounts)) - 0.5) * 0.1 + (i-1)*0.15;
                scatter(xPos, individualCounts, 50, color, 'o', 'filled', 'MarkerFaceAlpha', 0.3);
                
                % Store mean and SEM
                genotypeMeans_phase(d) = dayStats.(genField).(dayKey).(phase).mean_count;
                genotypeSEMs_phase(d) = dayStats.(genField).(dayKey).(phase).sem_count;
            end
        end
        
        % Plot means with error bars on top
        h_phase = errorbar(dayNumbers + (i-1)*0.15, genotypeMeans_phase, genotypeSEMs_phase, 'Color', color, 'LineStyle', 'none', ...
            'Marker', 'o', 'MarkerFaceColor', color, 'MarkerSize', 10, 'LineWidth', 2, 'CapSize', 6);
        
        % Add to legend arrays
        plotHandles_phase = [plotHandles_phase, h_phase];
        legendTexts_phase{end+1} = sprintf('%s (Mean ± SEM)', genotypeName);
    end
    
    title([phaseLabel ' Sleep Bout Counts Across Experimental Days']);
    xlabel('Experimental Day');
    ylabel('Number of Sleep Bouts');
    legend(plotHandles_phase, legendTexts_phase, 'Location', 'best');
    grid on;
    hold off;
    
    % Save figure
    saveas(gcf, fullfile(outputFolder, [strrep(phaseLabel, ' ', '_') '_Sleep_Bout_Counts_Across_Days.fig']));
    saveas(gcf, fullfile(outputFolder, [strrep(phaseLabel, ' ', '_') '_Sleep_Bout_Counts_Across_Days.png']));
end

fprintf('Per-day analysis complete!\n');

%% ZT hour analysis with dynamic genotypes
fprintf('\n==========================================\n');
fprintf('Beginning Zeitgeber Time (ZT) Analysis\n');
fprintf('==========================================\n');

% Define ZT hours (ZT0 = 6:00, ZT12 = 18:00)
numZTHours = 24;
ztLabels = cell(numZTHours, 1);
for i = 1:numZTHours
    ztLabels{i} = sprintf('ZT%d', i-1);
end

% Create structures to store data by ZT hour for each genotype and mouse
ztData = struct();
ztDataByDay = struct();

% Initialize ZT data structures for each genotype
for i = 1:numGenotypes
    genField = genotypeFieldNames{i};
    ztData.(genField) = struct('counts', zeros(0, numZTHours), 'totalDurations', zeros(0, numZTHours), 'avgDurations', zeros(0, numZTHours));
    ztDataByDay.(genField) = struct();
end

% Process each CSV file for ZT hour analysis
fprintf('Processing files for ZT hour analysis...\n');

for i = 1:length(csvFiles)
    csvPath = fullfile(csvFiles(i).folder, csvFiles(i).name);
    currentMouseID = mouseIDs{i};
    
    % Determine genotype based on user input
    mouseGenotype = '';
    genotypeDisplayName = '';
    for j = 1:numGenotypes
        genField = genotypeFieldNames{j};
        if ismember(currentMouseID, genotypeAssignments.(genField))
            mouseGenotype = genField;
            genotypeDisplayName = genotypeNames{j};
            break;
        end
    end
    
    if isempty(mouseGenotype)
        continue;
    end
    
    fprintf('Processing ZT data for mouse %s (genotype: %s)...\n', currentMouseID, genotypeDisplayName);
    
    try
        % Read the CSV file with all necessary columns
        opts = detectImportOptions(csvPath);
        
        % Debug: Print column names to verify what's in the file
        fprintf('  CSV columns: %s\n', strjoin(opts.VariableNames, ', '));
        
        % Make sure we have at least 3 columns
        if size(opts.VariableNames, 2) >= 3
            % Select timestamp (first column) and sleep bout duration (third column)
            opts.SelectedVariableNames = opts.VariableNames([1, 3]);
            
            % Set variable types appropriately
            opts = setvartype(opts, opts.SelectedVariableNames{1}, 'string'); % timestamp as string
            opts = setvartype(opts, opts.SelectedVariableNames{2}, 'double'); % duration as double
            
            % Read the data
            sleepData = readtable(csvPath, opts);
            
            if ~isempty(sleepData)
                % Extract sleep bout durations
                sleepDurations = sleepData.(opts.SelectedVariableNames{2});
                % Remove NaN values
                validIdx = ~isnan(sleepDurations);
                sleepDurations = sleepDurations(validIdx);
                
                % Extract timestamps for valid durations
                timestamps = sleepData.(opts.SelectedVariableNames{1})(validIdx);
                
                fprintf('  Found %d valid sleep bouts with timestamps\n', length(sleepDurations));
                % Debug output for the first few timestamps
                if ~isempty(timestamps)
                    fprintf('  First few timestamps: ');
                    for t = 1:min(3, length(timestamps))
                        fprintf('%s, ', string(timestamps(t)));
                    end
                    fprintf('\n');
                end
                
                % Group data by ZT hour
                ztCounts = zeros(1, numZTHours);
                ztTotalDurations = zeros(1, numZTHours);
                
                % Group data by ZT hour AND day
                ztCountsByDay = zeros(numDays, numZTHours);
                ztTotalDurationsByDay = zeros(numDays, numZTHours);
                
                % Process each timestamp and assign to ZT hour
                processedTimestamps = 0;
                failedTimestamps = 0;
                
                for t = 1:length(timestamps)
                    try
                        % Parse timestamp format (e.g., "9/19/23 6:01")
                        timestamp = timestamps(t);
                        
                        % Extract hour
                        timeParts = split(timestamp, ' ');
                        if length(timeParts) >= 2
                            timeStr = timeParts{2};
                            hourParts = split(timeStr, ':');
                            hour = str2double(hourParts{1});
                            
                            % Convert hour to ZT hour (ZT0 = 6:00, ZT12 = 18:00)
                            ztHour = mod(hour - 6, 24) + 1; % +1 for MATLAB 1-based indexing
                            
                            % Determine which experimental day this timestamp belongs to
                            timestampStr = char(timestamp);
                            dateParts = split(timestampStr, ' ');
                            datePart = dateParts{1}; % "2023-09-19"
                            timePart = dateParts{2}; % "06:07:00"
                            
                            % Extract hour from time part
                            timePartsDetailed = split(timePart, ':');
                            hourNum = str2double(timePartsDetailed{1});
                            
                            % Convert to datenum for easier manipulation
                            fullDatenum = datenum(timestampStr, 'yyyy-mm-dd HH:MM:SS');
                            
                            % Calculate experimental day using datenum arithmetic
                            if hourNum < 6
                                % Before 6 AM - belongs to previous day's experiment
                                expDayNum = floor(fullDatenum) - 1;
                            else
                                % 6 AM or later - belongs to current day's experiment
                                expDayNum = floor(fullDatenum);
                            end
                            
                            % Find which day this corresponds to
                            dayIdx = find(allDays == expDayNum);
                            if ~isempty(dayIdx)
                                dayNumber = dayIdx(1);
                                
                                % Increment count and add duration for this ZT hour and day
                                ztCountsByDay(dayNumber, ztHour) = ztCountsByDay(dayNumber, ztHour) + 1;
                                ztTotalDurationsByDay(dayNumber, ztHour) = ztTotalDurationsByDay(dayNumber, ztHour) + sleepDurations(t);
                            end
                            
                            % Increment overall count and add duration for this ZT hour
                            ztCounts(ztHour) = ztCounts(ztHour) + 1;
                            ztTotalDurations(ztHour) = ztTotalDurations(ztHour) + sleepDurations(t);
                            processedTimestamps = processedTimestamps + 1;
                        else
                            fprintf('    Warning: Unable to parse timestamp: %s\n', timestamp);
                            failedTimestamps = failedTimestamps + 1;
                        end
                    catch e
                        fprintf('    Warning: Error processing timestamp %s for ZT analysis: %s\n', timestamp, e.message);
                        failedTimestamps = failedTimestamps + 1;
                    end
                end
                
                fprintf('  Successfully processed %d of %d timestamps (%d failed)\n', processedTimestamps, length(timestamps), failedTimestamps);
                
                % Calculate average bout duration for each ZT hour
                ztAvgDurations = zeros(1, numZTHours);
                for h = 1:numZTHours
                    if ztCounts(h) > 0
                        ztAvgDurations(h) = ztTotalDurations(h) / ztCounts(h);
                    else
                        ztAvgDurations(h) = 0;
                    end
                end
                
                % Store the ZT data for this mouse
                prevSize = size(ztData.(mouseGenotype).counts, 1);
                ztData.(mouseGenotype).counts(prevSize+1,:) = ztCounts;
                ztData.(mouseGenotype).totalDurations(prevSize+1,:) = ztTotalDurations;
                ztData.(mouseGenotype).avgDurations(prevSize+1,:) = ztAvgDurations;
                
                % Store the per-day ZT data for this mouse
                mouseFieldName = makeValidFieldName(currentMouseID);
                ztDataByDay.(mouseGenotype).(mouseFieldName) = struct();
                ztDataByDay.(mouseGenotype).(mouseFieldName).counts = ztCountsByDay;
                ztDataByDay.(mouseGenotype).(mouseFieldName).totalDurations = ztTotalDurationsByDay;
                
                % Calculate average durations by day
                ztAvgDurationsByDay = zeros(numDays, numZTHours);
                for d = 1:numDays
                    for h = 1:numZTHours
                        if ztCountsByDay(d, h) > 0
                            ztAvgDurationsByDay(d, h) = ztTotalDurationsByDay(d, h) / ztCountsByDay(d, h);
                        end
                    end
                end
                ztDataByDay.(mouseGenotype).(mouseFieldName).avgDurations = ztAvgDurationsByDay;
                
                fprintf('  Added ZT data for mouse %s\n', currentMouseID);
                fprintf('  ZT summary for mouse %s:\n', currentMouseID);
                fprintf('    Total sleep bouts: %d\n', sum(ztCounts));
                fprintf('    Total sleep duration: %.1f seconds\n', sum(ztTotalDurations));
            else
                fprintf('  Warning: No sleep bout data found in file for mouse %s\n', currentMouseID);
            end
        else
            fprintf('  Warning: CSV file for mouse %s does not have at least 3 columns\n', currentMouseID);
        end
    catch e
        fprintf('  Error processing ZT data for mouse %s: %s\n', currentMouseID, e.message);
    end
end

%% Calculate ZT statistics by genotype
% Initialize structures for mean and std values
ztMean = struct();
ztStd = struct();
ztSem = struct(); % Standard error of the mean

for i = 1:numGenotypes
    genField = genotypeFieldNames{i};
    genotypeName = genotypeNames{i};
    
    if isfield(ztData, genField)
        % Calculate statistics for bout counts
        if ~isempty(ztData.(genField).counts)
            ztMean.(genField).counts = mean(ztData.(genField).counts, 1);
            ztStd.(genField).counts = std(ztData.(genField).counts, 0, 1);
            ztSem.(genField).counts = ztStd.(genField).counts / sqrt(size(ztData.(genField).counts, 1));
        else
            fprintf('No ZT bout count data found for genotype: %s\n', genotypeName);
            ztMean.(genField).counts = zeros(1, numZTHours);
            ztStd.(genField).counts = zeros(1, numZTHours);
            ztSem.(genField).counts = zeros(1, numZTHours);
        end
        
        % Calculate statistics for total durations
        if ~isempty(ztData.(genField).totalDurations)
            ztMean.(genField).totalDurations = mean(ztData.(genField).totalDurations, 1);
            ztStd.(genField).totalDurations = std(ztData.(genField).totalDurations, 0, 1);
            ztSem.(genField).totalDurations = ztStd.(genField).totalDurations / sqrt(size(ztData.(genField).totalDurations, 1));
        else
            fprintf('No ZT total duration data found for genotype: %s\n', genotypeName);
            ztMean.(genField).totalDurations = zeros(1, numZTHours);
            ztStd.(genField).totalDurations = zeros(1, numZTHours);
            ztSem.(genField).totalDurations = zeros(1, numZTHours);
        end
        
        % Calculate statistics for average durations
        if ~isempty(ztData.(genField).avgDurations)
            ztMean.(genField).avgDurations = mean(ztData.(genField).avgDurations, 1);
            ztStd.(genField).avgDurations = std(ztData.(genField).avgDurations, 0, 1);
            ztSem.(genField).avgDurations = ztStd.(genField).avgDurations / sqrt(size(ztData.(genField).avgDurations, 1));
        else
            fprintf('No ZT average duration data found for genotype: %s\n', genotypeName);
            ztMean.(genField).avgDurations = zeros(1, numZTHours);
            ztStd.(genField).avgDurations = zeros(1, numZTHours);
            ztSem.(genField).avgDurations = zeros(1, numZTHours);
        end
    end
end

%% Calculate per-day ZT statistics by genotype
fprintf('Calculating per-day ZT statistics...\n');

% Initialize structures for per-day ZT statistics
ztStatsByDay = struct();

for i = 1:numGenotypes
    genField = genotypeFieldNames{i};
    ztStatsByDay.(genField) = struct();
end

for d = 1:numDays
    dayKey = sprintf('day%d', d);
    
    for i = 1:numGenotypes
        genField = genotypeFieldNames{i};
        
        % Initialize storage for this day and genotype
        ztStatsByDay.(genField).(dayKey) = struct();
        ztStatsByDay.(genField).(dayKey).counts = struct('data', [], 'mean', zeros(1, numZTHours), 'std', zeros(1, numZTHours), 'sem', zeros(1, numZTHours));
        ztStatsByDay.(genField).(dayKey).totalDurations = struct('data', [], 'mean', zeros(1, numZTHours), 'std', zeros(1, numZTHours), 'sem', zeros(1, numZTHours));
        ztStatsByDay.(genField).(dayKey).avgDurations = struct('data', [], 'mean', zeros(1, numZTHours), 'std', zeros(1, numZTHours), 'sem', zeros(1, numZTHours));
        
        % Collect data from all mice of this genotype for this day
        if isfield(ztDataByDay, genField)
            mouseFields = fieldnames(ztDataByDay.(genField));
            
            for m = 1:length(mouseFields)
                mouseFieldName = mouseFields{m};
                mouseData = ztDataByDay.(genField).(mouseFieldName);
                
                % Add this mouse's data for this day
                ztStatsByDay.(genField).(dayKey).counts.data = [ztStatsByDay.(genField).(dayKey).counts.data; mouseData.counts(d, :)];
                ztStatsByDay.(genField).(dayKey).totalDurations.data = [ztStatsByDay.(genField).(dayKey).totalDurations.data; mouseData.totalDurations(d, :)];
                ztStatsByDay.(genField).(dayKey).avgDurations.data = [ztStatsByDay.(genField).(dayKey).avgDurations.data; mouseData.avgDurations(d, :)];
            end
            
            % Calculate statistics if we have data
            if ~isempty(ztStatsByDay.(genField).(dayKey).counts.data)
                ztStatsByDay.(genField).(dayKey).counts.mean = mean(ztStatsByDay.(genField).(dayKey).counts.data, 1);
                ztStatsByDay.(genField).(dayKey).counts.std = std(ztStatsByDay.(genField).(dayKey).counts.data, 0, 1);
                ztStatsByDay.(genField).(dayKey).counts.sem = ztStatsByDay.(genField).(dayKey).counts.std / sqrt(size(ztStatsByDay.(genField).(dayKey).counts.data, 1));
                
                ztStatsByDay.(genField).(dayKey).totalDurations.mean = mean(ztStatsByDay.(genField).(dayKey).totalDurations.data, 1);
                ztStatsByDay.(genField).(dayKey).totalDurations.std = std(ztStatsByDay.(genField).(dayKey).totalDurations.data, 0, 1);
                ztStatsByDay.(genField).(dayKey).totalDurations.sem = ztStatsByDay.(genField).(dayKey).totalDurations.std / sqrt(size(ztStatsByDay.(genField).(dayKey).totalDurations.data, 1));
                
                ztStatsByDay.(genField).(dayKey).avgDurations.mean = mean(ztStatsByDay.(genField).(dayKey).avgDurations.data, 1);
                ztStatsByDay.(genField).(dayKey).avgDurations.std = std(ztStatsByDay.(genField).(dayKey).avgDurations.data, 0, 1);
                ztStatsByDay.(genField).(dayKey).avgDurations.sem = ztStatsByDay.(genField).(dayKey).avgDurations.std / sqrt(size(ztStatsByDay.(genField).(dayKey).avgDurations.data, 1));
            end
        end
    end
end

%% Add ZT data to Excel export
fprintf('Adding ZT data to Excel export...\n');

% Create tables for ZT hour data
ztCountsTable = array2table(zeros(numZTHours, 0));
ztCountsTable.ZT_Hour = ztLabels;

ztTotalDurationsTable = array2table(zeros(numZTHours, 0));
ztTotalDurationsTable.ZT_Hour = ztLabels;

ztAvgDurationsTable = array2table(zeros(numZTHours, 0));
ztAvgDurationsTable.ZT_Hour = ztLabels;

% Add individual mouse data to tables
for i = 1:numGenotypes
    genField = genotypeFieldNames{i};
    
    if isfield(animalIDs, genField) && ~isempty(animalIDs.(genField))
        for m = 1:min(length(animalIDs.(genField)), size(ztData.(genField).counts, 1))
            mouseID = animalIDs.(genField){m};
            
            % Add counts data
            ztCountsTable.([mouseID '_Counts']) = ztData.(genField).counts(m,:)';
            
            % Add total durations data
            ztTotalDurationsTable.([mouseID '_Total']) = ztData.(genField).totalDurations(m,:)';
            
            % Add average durations data
            ztAvgDurationsTable.([mouseID '_Avg']) = ztData.(genField).avgDurations(m,:)';
        end
    end
end

% Add genotype statistics to tables
for i = 1:numGenotypes
    genField = genotypeFieldNames{i};
    display_name = genotypeNames{i}; % Use original genotype name
    
    if isfield(ztMean, genField)
        % Add bout counts statistics
        ztCountsTable.([display_name '_Mean']) = ztMean.(genField).counts';
        ztCountsTable.([display_name '_SD']) = ztStd.(genField).counts';
        ztCountsTable.([display_name '_SEM']) = ztSem.(genField).counts';
        
        % Add total durations statistics
        ztTotalDurationsTable.([display_name '_Mean']) = ztMean.(genField).totalDurations';
        ztTotalDurationsTable.([display_name '_SD']) = ztStd.(genField).totalDurations';
        ztTotalDurationsTable.([display_name '_SEM']) = ztSem.(genField).totalDurations';
        
        % Add average durations statistics
        ztAvgDurationsTable.([display_name '_Mean']) = ztMean.(genField).avgDurations';
        ztAvgDurationsTable.([display_name '_SD']) = ztStd.(genField).avgDurations';
        ztAvgDurationsTable.([display_name '_SEM']) = ztSem.(genField).avgDurations';
    end
end

% Write ZT data to Excel sheets
writetable(ztCountsTable, xlsFilePath, 'Sheet', 'ZT_Counts');
writetable(ztTotalDurationsTable, xlsFilePath, 'Sheet', 'ZT_Total_Duration');
writetable(ztAvgDurationsTable, xlsFilePath, 'Sheet', 'ZT_Avg_Duration');

fprintf('ZT data added to Excel file successfully.\n');

%% Create and write per-day ZT tables
fprintf('Creating per-day ZT data sheets...\n');

for d = 1:numDays
    dayKey = sprintf('day%d', d);
    dayDate = datestr(allDays(d), 'yyyy-mm-dd');
    
    fprintf('  Creating ZT data for Day %d (%s)...\n', d, dayDate);
    
    % Create tables for this day's ZT data
    dayZTCountsTable = array2table(zeros(numZTHours, 0));
    dayZTCountsTable.ZT_Hour = ztLabels;
    
    dayZTTotalDurTable = array2table(zeros(numZTHours, 0));
    dayZTTotalDurTable.ZT_Hour = ztLabels;
    
    % Add individual mouse data for this day
    allMouseIDs_ZT = {};
    for i = 1:numGenotypes
        genField = genotypeFieldNames{i};
        if isfield(animalIDs, genField)
            allMouseIDs_ZT = [allMouseIDs_ZT, animalIDs.(genField)];
        end
    end
    
    for m = 1:length(allMouseIDs_ZT)
        mouseID = allMouseIDs_ZT{m};
        mouseFieldName = makeValidFieldName(mouseID);
        
        % Find which genotype this mouse belongs to
        mouseGenotype = '';
        for i = 1:numGenotypes
            genField = genotypeFieldNames{i};
            if isfield(animalIDs, genField) && any(strcmp(animalIDs.(genField), mouseID))
                mouseGenotype = genField;
                break;
            end
        end
        
        if ~isempty(mouseGenotype)
            % Check if this mouse has ZT data
            if isfield(ztDataByDay, mouseGenotype) && isfield(ztDataByDay.(mouseGenotype), mouseFieldName)
                mouseZTData = ztDataByDay.(mouseGenotype).(mouseFieldName);
                
                % Add columns for this mouse
                dayZTCountsTable.([mouseID '_Counts']) = mouseZTData.counts(d, :)';
                dayZTTotalDurTable.([mouseID '_Total']) = mouseZTData.totalDurations(d, :)';
                
                fprintf('    Added ZT data for mouse %s\n', mouseID);
            else
                % Add zero columns if mouse has no ZT data
                zeroData = zeros(numZTHours, 1);
                dayZTCountsTable.([mouseID '_Counts']) = zeroData;
                dayZTTotalDurTable.([mouseID '_Total']) = zeroData;
            end
        end
    end
    
    % Add genotype statistics for this day
    for i = 1:numGenotypes
        genField = genotypeFieldNames{i};
        display_name = genotypeNames{i};
        
        if isfield(ztStatsByDay, genField) && isfield(ztStatsByDay.(genField), dayKey)
            dayStats = ztStatsByDay.(genField).(dayKey);
            
            % Add counts statistics
            dayZTCountsTable.([display_name '_Mean']) = dayStats.counts.mean';
            dayZTCountsTable.([display_name '_SD']) = dayStats.counts.std';
            dayZTCountsTable.([display_name '_SEM']) = dayStats.counts.sem';
            
            % Add total durations statistics
            dayZTTotalDurTable.([display_name '_Mean']) = dayStats.totalDurations.mean';
            dayZTTotalDurTable.([display_name '_SD']) = dayStats.totalDurations.std';
            dayZTTotalDurTable.([display_name '_SEM']) = dayStats.totalDurations.sem';
        end
    end
    
    % Write day ZT tables to Excel
    sheetNameCounts = sprintf('ZT_Day_%d_Counts', d);
    sheetNameTotalDur = sprintf('ZT_Day_%d_Total_Dur', d);
    
    try
        writetable(dayZTCountsTable, xlsFilePath, 'Sheet', sheetNameCounts);
        writetable(dayZTTotalDurTable, xlsFilePath, 'Sheet', sheetNameTotalDur);
        
        fprintf('  Added ZT sheets for Day %d: %s, %s\n', d, sheetNameCounts, sheetNameTotalDur);
    catch writeErr
        fprintf('  Warning: Could not write ZT Day %d data: %s\n', d, writeErr.message);
    end
end

%% Create ZT hour plots with dynamic genotypes
fprintf('Generating ZT hour plots...\n');

% Define plot types
ztPlotTypes = {'counts', 'totalDurations', 'avgDurations'};
ztPlotLabels = {'Sleep Bout Counts', 'Total Sleep Duration', 'Average Bout Duration'};
ztYLabels = {'Number of Bouts', 'Total Duration (seconds)', 'Average Duration (seconds)'};

% Create plots for each ZT data type
for p = 1:length(ztPlotTypes)
    plotType = ztPlotTypes{p};
    plotLabel = ztPlotLabels{p};
    yLabel = ztYLabels{p};
    
    %% Dot Plot (Individual mice with mean)
    figure('Name', ['ZT ' plotLabel ' - Dot Plot'], 'Position', [100, 100, 800, 500]);
    hold on;
    
    % Add light/dark background shading
    lightPhaseColor = [0.9, 0.9, 0.8];
    darkPhaseColor = [0.8, 0.8, 0.9];
    
    % Set axis limits first so background covers correctly
    xlim([-0.5, 23.5]);
    ylim auto;
    yl = ylim;
    
    % Draw background rectangles for light/dark phases
    rectangle('Position', [0, yl(1), 12, yl(2)-yl(1)], ...
        'FaceColor', lightPhaseColor, 'EdgeColor', 'none');
    rectangle('Position', [12, yl(1), 12, yl(2)-yl(1)], ...
        'FaceColor', darkPhaseColor, 'EdgeColor', 'none');
    
    % Initialize handles for legend
    plotHandles = [];
    legendTexts = {};
    
    % Plot data for each genotype
    for i = 1:numGenotypes
        genField = genotypeFieldNames{i};
        genotypeName = genotypeNames{i};
        color = genotypeColors(i, :);
        
        % Plot individual data points
        if isfield(ztData, genField) && ~isempty(ztData.(genField).(plotType))
            for j = 1:size(ztData.(genField).(plotType), 1)
                scatter((0:23) + (i-1)*0.1, ztData.(genField).(plotType)(j,:), 50, color, 'o', 'filled', 'MarkerFaceAlpha', 0.3);
            end
            
            % Plot mean with error bars and capture handle
            h = errorbar((0:23) + (i-1)*0.1, ztMean.(genField).(plotType), ztSem.(genField).(plotType), ...
                'Color', color, 'LineWidth', 2, 'Marker', 'o', 'MarkerFaceColor', color, 'MarkerSize', 10);
            
            % Add to legend arrays
            plotHandles = [plotHandles, h];
            legendTexts{end+1} = sprintf('%s (Mean ± SEM)', genotypeName);
        end
    end
    
    % Bring the lines to front
    uistack(findobj(gca, 'Type', 'line'), 'top');
    
    % Add labels and legend
    title(['ZT Hour ' plotLabel]);
    xlabel('ZT Hour (ZT0 = 6:00, ZT12 = 18:00)');
    ylabel(yLabel);
    set(gca, 'XTick', 0:2:23);
    
    % Create legend with handles
    if ~isempty(plotHandles)
        legend(plotHandles, legendTexts, 'Location', 'best');
    end
    
    grid on;
    hold off;
    
    % Save the figure
    saveas(gcf, fullfile(outputFolder, ['ZT_' strrep(plotLabel, ' ', '_') '_Dot.fig']));
    saveas(gcf, fullfile(outputFolder, ['ZT_' strrep(plotLabel, ' ', '_') '_Dot.png']));
end

fprintf('ZT hour analysis complete!\n');

%% Create ZT plots across experimental days
fprintf('Creating ZT bout counts across experimental days plots...\n');

% Define colors for each day
if numDays <= 5
    dayColors = [
        0.2, 0.4, 0.8;  % Day 1 - Blue
        0.8, 0.2, 0.2;  % Day 2 - Red  
        0.2, 0.8, 0.2;  % Day 3 - Green
        0.8, 0.4, 0.2;  % Day 4 - Orange
        0.6, 0.2, 0.8;  % Day 5 - Purple
    ];
    dayColors = dayColors(1:numDays, :);
else
    % Use lines colormap for many days
    dayColors = colormap(lines(numDays));
end

% 1. ZT bout counts across days by genotype
for i = 1:numGenotypes
    genField = genotypeFieldNames{i};
    display_name = genotypeNames{i};
    
    figure('Name', ['ZT Bout Counts Across Days - ' display_name], 'Position', [100, 100, 1200, 800]);
    hold on;
    
    % Add light/dark background shading
    lightPhaseColor = [0.9, 0.9, 0.8];
    darkPhaseColor = [0.8, 0.8, 0.9];
    
    % Set axis limits and draw background
    xlim([-0.5, 23.5]);
    ylim auto;
    yl = ylim;
    
    rectangle('Position', [0, yl(1), 12, yl(2)-yl(1)], ...
        'FaceColor', lightPhaseColor, 'EdgeColor', 'none');
    rectangle('Position', [12, yl(1), 12, yl(2)-yl(1)], ...
        'FaceColor', darkPhaseColor, 'EdgeColor', 'none');
    
    % Plot data for each day
    plotHandles = [];
    legendTexts = {};
    
    for d = 1:numDays
        dayKey = sprintf('day%d', d);
        
        if isfield(ztStatsByDay, genField) && isfield(ztStatsByDay.(genField), dayKey) && ...
           ~isempty(ztStatsByDay.(genField).(dayKey).counts.mean)
            
            meanData = ztStatsByDay.(genField).(dayKey).counts.mean;
            semData = ztStatsByDay.(genField).(dayKey).counts.sem;
            
            h = errorbar(0:23, meanData, semData, ...
                'Color', dayColors(d, :), 'LineWidth', 2, 'Marker', 'o', ...
                'MarkerFaceColor', dayColors(d, :), 'MarkerSize', 6);
            
            plotHandles = [plotHandles, h];
            legendTexts{end+1} = sprintf('Day %d', d);
        end
    end
    
    % Bring plot lines to front
    uistack(findobj(gca, 'Type', 'line'), 'top');
    
    title(['ZT Bout Counts Across Days - ' display_name ' Mice']);
    xlabel('ZT Hour (ZT0 = 6:00, ZT12 = 18:00)');
    ylabel('Number of Sleep Bouts (Mean ± SEM)');
    set(gca, 'XTick', 0:2:23);
    
    if ~isempty(plotHandles)
        legend(plotHandles, legendTexts, 'Location', 'best');
    end
    
    grid on;
    hold off;
    
    % Save figure
    saveas(gcf, fullfile(outputFolder, ['ZT_Counts_Across_Days_' strrep(display_name, '-', '_') '.fig']));
    saveas(gcf, fullfile(outputFolder, ['ZT_Counts_Across_Days_' strrep(display_name, '-', '_') '.png']));
end

% 2. Comparison plot: Wild-type vs Mutant for each day
for d = 1:numDays
    dayKey = sprintf('day%d', d);
    dayDate = datestr(allDays(d), 'yyyy-mm-dd');
    
    figure('Name', ['ZT Bout Counts Day ' num2str(d) ' - Genotype Comparison'], 'Position', [100, 100, 1000, 600]);
    hold on;
    
    % Add light/dark background shading
    lightPhaseColor = [0.9, 0.9, 0.8];
    darkPhaseColor = [0.8, 0.8, 0.9];
    
    xlim([-0.5, 23.5]);
    ylim auto;
    yl = ylim;
    
    rectangle('Position', [0, yl(1), 12, yl(2)-yl(1)], ...
        'FaceColor', lightPhaseColor, 'EdgeColor', 'none');
    rectangle('Position', [12, yl(1), 12, yl(2)-yl(1)], ...
        'FaceColor', darkPhaseColor, 'EdgeColor', 'none');
    
    % Plot all genotypes
    plotHandles = [];
    legendTexts = {};
    
    for i = 1:numGenotypes
        genField = genotypeFieldNames{i};
        display_name = genotypeNames{i};
        color = genotypeColors(i, :);
        
        if isfield(ztStatsByDay, genField) && isfield(ztStatsByDay.(genField), dayKey) && ...
           ~isempty(ztStatsByDay.(genField).(dayKey).counts.mean)
            
            meanData = ztStatsByDay.(genField).(dayKey).counts.mean;
            semData = ztStatsByDay.(genField).(dayKey).counts.sem;
            
            h = errorbar(0:23, meanData, semData, ...
                'Color', color, 'LineWidth', 2, 'Marker', 'o', ...
                'MarkerFaceColor', color, 'MarkerSize', 8);
            
            plotHandles = [plotHandles, h];
            legendTexts{end+1} = [display_name ' (Mean ± SEM)'];
        end
    end
    
    % Bring plot lines to front
    uistack(findobj(gca, 'Type', 'line'), 'top');
    
    title(['ZT Bout Counts - Day ' num2str(d) ' (' dayDate ') - Genotype Comparison']);
    xlabel('ZT Hour (ZT0 = 6:00, ZT12 = 18:00)');
    ylabel('Number of Sleep Bouts (Mean ± SEM)');
    set(gca, 'XTick', 0:2:23);
    
    if ~isempty(plotHandles)
        legend(plotHandles, legendTexts, 'Location', 'best');
    end
    
    grid on;
    hold off;
    
    % Save figure
    saveas(gcf, fullfile(outputFolder, ['ZT_Counts_Day_' num2str(d) '_Genotype_Comparison.fig']));
    saveas(gcf, fullfile(outputFolder, ['ZT_Counts_Day_' num2str(d) '_Genotype_Comparison.png']));
end

% 3. Individual mouse ZT patterns across days (create subplots for each mouse)
allMouseIDs_plot = {};
for i = 1:numGenotypes
    genField = genotypeFieldNames{i};
    if isfield(animalIDs, genField)
        allMouseIDs_plot = [allMouseIDs_plot, animalIDs.(genField)];
    end
end

% Create individual mouse plots (4 mice per figure)
micePerFig = 4;
numFigs = ceil(length(allMouseIDs_plot) / micePerFig);

for figIdx = 1:numFigs
    figure('Name', ['Individual Mouse ZT Patterns Across Days - Set ' num2str(figIdx)], 'Position', [100, 100, 1400, 1000]);
    
    startIdx = (figIdx - 1) * micePerFig + 1;
    endIdx = min(figIdx * micePerFig, length(allMouseIDs_plot));
    
    for plotIdx = 1:(endIdx - startIdx + 1)
        mouseIdx = startIdx + plotIdx - 1;
        mouseID = allMouseIDs_plot{mouseIdx};
        mouseFieldName = makeValidFieldName(mouseID);
        
        % Find which genotype this mouse belongs to
        mouseGenotype = '';
        mouseColor = [];
        for i = 1:numGenotypes
            genField = genotypeFieldNames{i};
            if isfield(animalIDs, genField) && any(strcmp(animalIDs.(genField), mouseID))
                mouseGenotype = genField;
                mouseColor = genotypeColors(i, :);
                break;
            end
        end
        
        if isempty(mouseGenotype)
            continue;
        end
        
        subplot(2, 2, plotIdx);
        hold on;
        
        % Add light/dark background shading
        lightPhaseColor = [0.9, 0.9, 0.8];
        darkPhaseColor = [0.8, 0.8, 0.9];
        
        xlim([-0.5, 23.5]);
        ylim auto;
        yl = ylim;
        
        rectangle('Position', [0, yl(1), 12, yl(2)-yl(1)], ...
            'FaceColor', lightPhaseColor, 'EdgeColor', 'none');
        rectangle('Position', [12, yl(1), 12, yl(2)-yl(1)], ...
            'FaceColor', darkPhaseColor, 'EdgeColor', 'none');
        
        % Plot this mouse's data across days
        if isfield(ztDataByDay, mouseGenotype) && isfield(ztDataByDay.(mouseGenotype), mouseFieldName)
            mouseZTData = ztDataByDay.(mouseGenotype).(mouseFieldName);
            
            for d = 1:numDays
                plot(0:23, mouseZTData.counts(d, :), ...
                    'Color', dayColors(d, :), 'LineWidth', 1.5, 'Marker', 'o', ...
                    'MarkerFaceColor', dayColors(d, :), 'MarkerSize', 4);
            end
        end
        
        % Bring plot lines to front
        uistack(findobj(gca, 'Type', 'line'), 'top');
        
        % Find genotype display name
        genotypeDisplayName = '';
        for i = 1:numGenotypes
            genField = genotypeFieldNames{i};
            if strcmp(genField, mouseGenotype)
                genotypeDisplayName = genotypeNames{i};
                break;
            end
        end
        
        title(['Mouse ' mouseID ' (' genotypeDisplayName ')']);
        xlabel('ZT Hour');
        ylabel('Sleep Bout Counts');
        set(gca, 'XTick', 0:4:23);
        grid on;
        hold off;
        
        % Add legend only to first subplot
        if plotIdx == 1
            legendTexts_mouse = {};
            for d = 1:numDays
                legendTexts_mouse{end+1} = sprintf('Day %d', d);
            end
            legend(legendTexts_mouse, 'Location', 'best', 'FontSize', 8);
        end
    end
    
    % Save figure
    saveas(gcf, fullfile(outputFolder, ['Individual_Mouse_ZT_Patterns_Set_' num2str(figIdx) '.fig']));
    saveas(gcf, fullfile(outputFolder, ['Individual_Mouse_ZT_Patterns_Set_' num2str(figIdx) '.png']));
end

fprintf('ZT across days analysis complete!\n');

% Restore normal figure visibility
set(0, 'DefaultFigureVisible', 'on');

fprintf('ZT analysis complete!\n');
fprintf('All figures saved to %s\n', outputFolder);
fprintf('Analysis complete!\n');