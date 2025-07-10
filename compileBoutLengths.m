% This script:
% 1. Finds all CSV files in input folder (no subfolders)
% 2. Groups mice by genotype (wild-type or mutant)
% 3. Retrieves epoch durations for each mouse from CSV files
% 4. Sorts sleep bouts into light phase (6:00-18:00) and dark phase (18:00-6:00)
% 5. Bins the epochs by their duration (2s, 4s, 8s, 16s, 32s, 64s, 128s, 256s, 512s)
% 6. Performs ZT (Zeitgeber Time) hour analysis, where ZT0 = 6:00, ZT12 = 18:00
% 7. SEPARATES DATA BY EXPERIMENTAL DAY (Day 1, Day 2, etc.) - days run from 6AM to 6AM
% 8. Computes averages and standard deviations by genotype (overall, light phase, dark phase, ZT hour, per day)
% 9. Plots results as dot graphs with error bars for all data categories

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
genotypes = {'wild-type', 'mutant'}; % Define genotypes 

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

% Create structures to store binned epoch counts by genotype
epochCounts = struct();
epochCounts.wild_type = struct('sleep', zeros(0, numBins), 'light', zeros(0, numBins), 'dark', zeros(0, numBins));
epochCounts.mutant = struct('sleep', zeros(0, numBins), 'light', zeros(0, numBins), 'dark', zeros(0, numBins));

% Create structure to store animal IDs by genotype
animalIDs = struct();
animalIDs.wild_type = {};
animalIDs.mutant = {};

%% Initialize per-day data structures
% These will store data separated by experimental day
dayData = struct();
dayData.wild_type = struct();
dayData.mutant = struct();

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

% Prompt user to identify wild-type mice
fprintf('Please enter the IDs of wild-type mice, separated by commas:\n');
wtInput = input('Wild-type mice: ', 's');
% Process the input - split by commas and trim whitespace
wtMice = strtrim(split(wtInput, ','));
for i = 1:length(wtMice)
    wtMice{i} = strtrim(wtMice{i});
end

% Prompt user to identify mutant mice
fprintf('Please enter the IDs of mutant mice, separated by commas:\n');
mutInput = input('Mutant mice: ', 's');
% Process the input - split by commas and trim whitespace
mutMice = strtrim(split(mutInput, ','));
for i = 1:length(mutMice)
    mutMice{i} = strtrim(mutMice{i});
end

% Verify inputs and warn about any unlisted mice
allListedMice = [wtMice; mutMice];
unlistedMice = setdiff(mouseIDs, allListedMice);
if ~isempty(unlistedMice)
    fprintf('\nWarning: The following mice were not assigned to any genotype and will be skipped:\n');
    for i = 1:length(unlistedMice)
        fprintf('  %s\n', unlistedMice{i});
    end
    fprintf('\n');
end

% Make sure there's no overlap between wild-type and mutant mice
commonMice = intersect(wtMice, mutMice);
if ~isempty(commonMice)
    fprintf('\nERROR: The following mice were assigned to both genotypes:\n');
    for i = 1:length(commonMice)
        fprintf('  %s\n', commonMice{i});
    end
    error('Mice cannot be assigned to both genotypes. Please restart the script and provide non-overlapping lists.');
end

fprintf('Processing will begin with:\n');
fprintf('  Wild-type mice: %s\n', strjoin(wtMice, ', '));
fprintf('  Mutant mice: %s\n\n', strjoin(mutMice, ', '));

%% Process each CSV file
% Create a structure to store total epoch counts per mouse
totalCounts = struct();
totalCounts.wild_type = struct('sleep', [], 'light', [], 'dark', []);
totalCounts.mutant = struct('sleep', [], 'light', [], 'dark', []);

% Match files to mice and process them
for i = 1:length(csvFiles)
    csvPath = fullfile(csvFiles(i).folder, csvFiles(i).name);
    currentMouseID = mouseIDs{i};
    
    % Determine genotype based on user input
    if ismember(currentMouseID, wtMice)
        mouseGenotype = 'wild_type';
    elseif ismember(currentMouseID, mutMice)
        mouseGenotype = 'mutant';
    else
        % Skip mice not assigned to any genotype
        fprintf('Skipping mouse %s (no genotype assigned)...\n', currentMouseID);
        continue;
    end
    
    fprintf('Processing mouse %s (genotype: %s)...\n', currentMouseID, strrep(mouseGenotype, '_', '-'));
    
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
                if ~isfield(totalCounts.(mouseGenotype), 'mouseIDs')
                    totalCounts.(mouseGenotype).mouseIDs = {};
                end
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

% Create structures to store per-day statistics
dayStats = struct();
dayStats.wild_type = struct();
dayStats.mutant = struct();

numDays = length(allDays);

% For each genotype, calculate stats for each day
for genotype = {'wild_type', 'mutant'}
    gen = genotype{1};
    
    % Initialize per-day storage
    for d = 1:numDays
        dayKey = sprintf('day%d', d);
        dayStats.(gen).(dayKey) = struct();
        dayStats.(gen).(dayKey).sleep = struct('counts', [], 'binned', []);
        dayStats.(gen).(dayKey).light = struct('counts', [], 'binned', []);
        dayStats.(gen).(dayKey).dark = struct('counts', [], 'binned', []);
    end
    
    % Process each mouse of this genotype
    if isfield(dayData, gen)
        mouseFields = fieldnames(dayData.(gen));
        for m = 1:length(mouseFields)
            mouseFieldName = mouseFields{m};
            mouseData = dayData.(gen).(mouseFieldName);
            
            % Process each day for this mouse
            mouseDayFields = fieldnames(mouseData);
            for d = 1:length(mouseDayFields)
                dayKey = mouseDayFields{d};
                
                % Get data for this day
                sleepBouts = mouseData.(dayKey).sleep;
                lightBouts = mouseData.(dayKey).light;
                darkBouts = mouseData.(dayKey).dark;
                
                % Store counts
                dayStats.(gen).(dayKey).sleep.counts = [dayStats.(gen).(dayKey).sleep.counts; length(sleepBouts)];
                dayStats.(gen).(dayKey).light.counts = [dayStats.(gen).(dayKey).light.counts; length(lightBouts)];
                dayStats.(gen).(dayKey).dark.counts = [dayStats.(gen).(dayKey).dark.counts; length(darkBouts)];
                
                % Bin the data
                sleepBinned = histcounts(sleepBouts, binEdges);
                lightBinned = histcounts(lightBouts, binEdges);
                darkBinned = histcounts(darkBouts, binEdges);
                
                % Store binned data
                dayStats.(gen).(dayKey).sleep.binned = [dayStats.(gen).(dayKey).sleep.binned; sleepBinned];
                dayStats.(gen).(dayKey).light.binned = [dayStats.(gen).(dayKey).light.binned; lightBinned];
                dayStats.(gen).(dayKey).dark.binned = [dayStats.(gen).(dayKey).dark.binned; darkBinned];
            end
        end
    end
    
    % Calculate means and SEMs for each day
    for d = 1:numDays
        dayKey = sprintf('day%d', d);
        
        % Sleep bouts
        if ~isempty(dayStats.(gen).(dayKey).sleep.counts)
            dayStats.(gen).(dayKey).sleep.mean_count = mean(dayStats.(gen).(dayKey).sleep.counts);
            dayStats.(gen).(dayKey).sleep.sem_count = std(dayStats.(gen).(dayKey).sleep.counts) / sqrt(length(dayStats.(gen).(dayKey).sleep.counts));
            dayStats.(gen).(dayKey).sleep.mean_binned = mean(dayStats.(gen).(dayKey).sleep.binned, 1);
            dayStats.(gen).(dayKey).sleep.sem_binned = std(dayStats.(gen).(dayKey).sleep.binned, 0, 1) / sqrt(size(dayStats.(gen).(dayKey).sleep.binned, 1));
        end
        
        % Light phase
        if ~isempty(dayStats.(gen).(dayKey).light.counts)
            dayStats.(gen).(dayKey).light.mean_count = mean(dayStats.(gen).(dayKey).light.counts);
            dayStats.(gen).(dayKey).light.sem_count = std(dayStats.(gen).(dayKey).light.counts) / sqrt(length(dayStats.(gen).(dayKey).light.counts));
            dayStats.(gen).(dayKey).light.mean_binned = mean(dayStats.(gen).(dayKey).light.binned, 1);
            dayStats.(gen).(dayKey).light.sem_binned = std(dayStats.(gen).(dayKey).light.binned, 0, 1) / sqrt(size(dayStats.(gen).(dayKey).light.binned, 1));
        end
        
        % Dark phase
        if ~isempty(dayStats.(gen).(dayKey).dark.counts)
            dayStats.(gen).(dayKey).dark.mean_count = mean(dayStats.(gen).(dayKey).dark.counts);
            dayStats.(gen).(dayKey).dark.sem_count = std(dayStats.(gen).(dayKey).dark.counts) / sqrt(length(dayStats.(gen).(dayKey).dark.counts));
            dayStats.(gen).(dayKey).dark.mean_binned = mean(dayStats.(gen).(dayKey).dark.binned, 1);
            dayStats.(gen).(dayKey).dark.sem_binned = std(dayStats.(gen).(dayKey).dark.binned, 0, 1) / sqrt(size(dayStats.(gen).(dayKey).dark.binned, 1));
        end
    end
end

%% Calculate statistics by genotype
% Initialize structures for mean and std values
meanCounts = struct();
stdCounts = struct();
semCounts = struct(); % Standard error of the mean

for genotype = {'wild_type', 'mutant'}
    gen = genotype{1};
    if isfield(epochCounts, gen)
        % Calculate mean and std for all sleep bouts
        if ~isempty(epochCounts.(gen).sleep)
            meanCounts.(gen).sleep = mean(epochCounts.(gen).sleep, 1);
            stdCounts.(gen).sleep = std(epochCounts.(gen).sleep, 0, 1);
            semCounts.(gen).sleep = stdCounts.(gen).sleep / sqrt(size(epochCounts.(gen).sleep, 1));
        else
            fprintf('No sleep bout data found for genotype: %s\n', gen);
            meanCounts.(gen).sleep = zeros(1, numBins);
            stdCounts.(gen).sleep = zeros(1, numBins);
            semCounts.(gen).sleep = zeros(1, numBins);
        end
        
        % Calculate mean and std for light phase bouts
        if ~isempty(epochCounts.(gen).light)
            meanCounts.(gen).light = mean(epochCounts.(gen).light, 1);
            stdCounts.(gen).light = std(epochCounts.(gen).light, 0, 1);
            semCounts.(gen).light = stdCounts.(gen).light / sqrt(size(epochCounts.(gen).light, 1));
        else
            fprintf('No light phase data found for genotype: %s\n', gen);
            meanCounts.(gen).light = zeros(1, numBins);
            stdCounts.(gen).light = zeros(1, numBins);
            semCounts.(gen).light = zeros(1, numBins);
        end
        
        % Calculate mean and std for dark phase bouts
        if ~isempty(epochCounts.(gen).dark)
            meanCounts.(gen).dark = mean(epochCounts.(gen).dark, 1);
            stdCounts.(gen).dark = std(epochCounts.(gen).dark, 0, 1);
            semCounts.(gen).dark = stdCounts.(gen).dark / sqrt(size(epochCounts.(gen).dark, 1));
        else
            fprintf('No dark phase data found for genotype: %s\n', gen);
            meanCounts.(gen).dark = zeros(1, numBins);
            stdCounts.(gen).dark = zeros(1, numBins);
            semCounts.(gen).dark = zeros(1, numBins);
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
for genotype = {'wild_type', 'mutant'}
    gen = genotype{1};
    
    if isfield(animalIDs, gen) && ~isempty(animalIDs.(gen))
        for m = 1:length(animalIDs.(gen))
            mouseID = animalIDs.(gen){m};
            fprintf('  Adding data for mouse: %s\n', mouseID);
            
            % Create data columns for all sleep phases
            if isfield(totalCounts, gen) && length(totalCounts.(gen).sleep) >= m
                sleepData = [totalCounts.(gen).sleep(m); epochCounts.(gen).sleep(m,:)'];
                lightData = [totalCounts.(gen).light(m); epochCounts.(gen).light(m,:)'];
                darkData = [totalCounts.(gen).dark(m); epochCounts.(gen).dark(m,:)'];
                
                % Add all columns for this mouse
                totalBinTable.([mouseID '_Sleep']) = sleepData;
                totalBinTable.([mouseID '_Light']) = lightData;
                totalBinTable.([mouseID '_Dark']) = darkData;
            end
        end
    end
end

% Add genotype statistics for all phases
for genotype = {'wild_type', 'mutant'}
    gen = genotype{1};
    display_name = strrep(gen, '_', '-'); % Convert 'wild_type' to 'wild-type'
    
    if isfield(totalCounts, gen)
        % Sleep statistics
        totalBinTable.([display_name '_Sleep_Mean']) = [mean(totalCounts.(gen).sleep); meanCounts.(gen).sleep'];
        totalBinTable.([display_name '_Sleep_SD']) = [std(totalCounts.(gen).sleep); stdCounts.(gen).sleep'];
        totalSEM_sleep = std(totalCounts.(gen).sleep)/sqrt(length(totalCounts.(gen).sleep));
        totalBinTable.([display_name '_Sleep_SEM']) = [totalSEM_sleep; semCounts.(gen).sleep'];
        
        % Light phase statistics
        totalBinTable.([display_name '_Light_Mean']) = [mean(totalCounts.(gen).light); meanCounts.(gen).light'];
        totalBinTable.([display_name '_Light_SD']) = [std(totalCounts.(gen).light); stdCounts.(gen).light'];
        totalSEM_light = std(totalCounts.(gen).light)/sqrt(length(totalCounts.(gen).light));
        totalBinTable.([display_name '_Light_SEM']) = [totalSEM_light; semCounts.(gen).light'];
        
        % Dark phase statistics
        totalBinTable.([display_name '_Dark_Mean']) = [mean(totalCounts.(gen).dark); meanCounts.(gen).dark'];
        totalBinTable.([display_name '_Dark_SD']) = [std(totalCounts.(gen).dark); stdCounts.(gen).dark'];
        totalSEM_dark = std(totalCounts.(gen).dark)/sqrt(length(totalCounts.(gen).dark));
        totalBinTable.([display_name '_Dark_SEM']) = [totalSEM_dark; semCounts.(gen).dark'];
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
    for genotype = {'wild_type', 'mutant'}
        gen = genotype{1};
        
        if isfield(animalIDs, gen) && ~isempty(animalIDs.(gen))
            for m = 1:length(animalIDs.(gen))
                mouseID = animalIDs.(gen){m};
                mouseFieldName = makeValidFieldName(mouseID);
                
                % Check if this mouse has data for this day
                if isfield(dayData, gen) && isfield(dayData.(gen), mouseFieldName)
                    mouseDayData = dayData.(gen).(mouseFieldName);
                    
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
    for genotype = {'wild_type', 'mutant'}
        gen = genotype{1};
        display_name = strrep(gen, '_', '-');
        
        if isfield(dayStats, gen) && isfield(dayStats.(gen), dayKey)
            % Add mean and SEM for sleep, light, and dark phases
            phases = {'sleep', 'light', 'dark'};
            
            for p = 1:length(phases)
                phase = phases{p};
                
                if isfield(dayStats.(gen).(dayKey), phase) && ...
                   isfield(dayStats.(gen).(dayKey).(phase), 'mean_count')
                    
                    % Total count + binned means
                    meanData = [dayStats.(gen).(dayKey).(phase).mean_count; ...
                               dayStats.(gen).(dayKey).(phase).mean_binned'];
                    
                    % Total SEM + binned SEMs  
                    semData = [dayStats.(gen).(dayKey).(phase).sem_count; ...
                              dayStats.(gen).(dayKey).(phase).sem_binned'];
                    
                    dayTable.([display_name '_' phase '_Mean']) = meanData;
                    dayTable.([display_name '_' phase '_SEM']) = semData;
                end
            end
        end
    end
    
    % Write day table to Excel with new sheet name
    sheetName = sprintf('Day_%d_Bin', d);  % Changed from 'Day_%d_GT_Bin' to 'Day_%d_Bin'
    xlsFilePath = fullfile(mainFolder, 'sleep_bout_analysis.xls');
    writetable(dayTable, xlsFilePath, 'Sheet', sheetName);
    
    fprintf('  Added Day %d data to Excel sheet: %s\n', d, sheetName);
end

% 2. Create a summary table mapping mice to genotypes
mouseGenotypeSummary = table();
mouseGenotypeSummary.MouseID = [animalIDs.wild_type(:); animalIDs.mutant(:)];
genotypes_list = [repmat({'wild-type'}, length(animalIDs.wild_type), 1); 
                  repmat({'mutant'}, length(animalIDs.mutant), 1)];
mouseGenotypeSummary.Genotype = genotypes_list;

% Write all tables to different sheets in the Excel file
xlsFilePath = fullfile(mainFolder, 'sleep_bout_analysis.xls');
fprintf('Writing combined sheets to Excel file: %s\n', xlsFilePath);

% Write sheets in desired order: Per-day sheets first, then combined sheets
% Note: Per-day sheets were already written in the loop above

% Write combined sheet - now includes all data (total + light/dark + genotype stats)
writetable(totalBinTable, xlsFilePath, 'Sheet', 'Total_Bin');
writetable(mouseGenotypeSummary, xlsFilePath, 'Sheet', 'Mouse_Genotype_Summary');

fprintf('Excel file created successfully.\n');

%% Create figures
% Turn off figure visibility to prevent windows from opening
set(0, 'DefaultFigureVisible', 'off');

% Colors for plotting
wtColor = [0, 0, 0.8]; % Blue for wild-type
mutColor = [0.8, 0, 0]; % Red for mutant

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
    
    % Plot individual data points for wild-type
    if isfield(epochCounts, 'wild_type') && ~isempty(epochCounts.wild_type.(category))
        for i = 1:size(epochCounts.wild_type.(category), 1)
            scatter(1:numBins, epochCounts.wild_type.(category)(i,:), 50, wtColor, 'o', 'filled', 'MarkerFaceAlpha', 0.3);
        end
        
        % Plot mean with error bars and capture handle
        h_wt = errorbar(1:numBins, meanCounts.wild_type.(category), semCounts.wild_type.(category), ...
            'Color', wtColor, 'LineWidth', 2, 'Marker', 'o', 'MarkerFaceColor', wtColor, 'MarkerSize', 10);
        
        % Add to legend arrays
        plotHandles = [plotHandles, h_wt];
        legendTexts{end+1} = 'Wild-type (Mean ± SEM)';
    end
    
    % Plot individual data points for mutant
    if isfield(epochCounts, 'mutant') && ~isempty(epochCounts.mutant.(category))
        for i = 1:size(epochCounts.mutant.(category), 1)
            scatter((1:numBins)+0.2, epochCounts.mutant.(category)(i,:), 50, mutColor, 'o', 'filled', 'MarkerFaceAlpha', 0.3);
        end
        
        % Plot mean with error bars and capture handle
        h_mut = errorbar((1:numBins)+0.2, meanCounts.mutant.(category), semCounts.mutant.(category), ...
            'Color', mutColor, 'LineWidth', 2, 'Marker', 'o', 'MarkerFaceColor', mutColor, 'MarkerSize', 10);
        
        % Add to legend arrays
        plotHandles = [plotHandles, h_mut];
        legendTexts{end+1} = 'Mutant (Mean ± SEM)';
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

%% Create per-day plots
fprintf('Creating per-day plots...\n');

% Plot sleep bout counts across days
figure('Name', 'Sleep Bout Counts Across Days', 'Position', [100, 100, 1000, 600]);
hold on;

dayNumbers = 1:numDays;
wtMeans = zeros(1, numDays);
wtSEMs = zeros(1, numDays);
mutMeans = zeros(1, numDays);
mutSEMs = zeros(1, numDays);

% Initialize handles for legend
plotHandles = [];
legendTexts = {};

% Plot individual data points for wild-type mice
for d = 1:numDays
    dayKey = sprintf('day%d', d);
    
    if isfield(dayStats.wild_type, dayKey) && ~isempty(dayStats.wild_type.(dayKey).sleep.counts)
        individualCounts = dayStats.wild_type.(dayKey).sleep.counts;
        % Add small random jitter to x-position for visibility
        xPos = d + (rand(size(individualCounts)) - 0.5) * 0.1;
        scatter(xPos, individualCounts, 50, wtColor, 'o', 'filled', 'MarkerFaceAlpha', 0.3);
        
        % Store mean and SEM
        wtMeans(d) = dayStats.wild_type.(dayKey).sleep.mean_count;
        wtSEMs(d) = dayStats.wild_type.(dayKey).sleep.sem_count;
    end
end

% Plot individual data points for mutant mice
for d = 1:numDays
    dayKey = sprintf('day%d', d);
    
    if isfield(dayStats.mutant, dayKey) && ~isempty(dayStats.mutant.(dayKey).sleep.counts)
        individualCounts = dayStats.mutant.(dayKey).sleep.counts;
        % Add small random jitter to x-position for visibility
        xPos = d + (rand(size(individualCounts)) - 0.5) * 0.1 + 0.2;
        scatter(xPos, individualCounts, 50, mutColor, 'o', 'filled', 'MarkerFaceAlpha', 0.3);
        
        % Store mean and SEM
        mutMeans(d) = dayStats.mutant.(dayKey).sleep.mean_count;
        mutSEMs(d) = dayStats.mutant.(dayKey).sleep.sem_count;
    end
end

% Plot means with error bars on top
h_wt = errorbar(dayNumbers, wtMeans, wtSEMs, 'Color', wtColor, 'LineStyle', 'none', ...
    'Marker', 'o', 'MarkerFaceColor', wtColor, 'MarkerSize', 10, 'LineWidth', 2, 'CapSize', 6);
h_mut = errorbar(dayNumbers + 0.2, mutMeans, mutSEMs, 'Color', mutColor, 'LineStyle', 'none', ...
    'Marker', 'o', 'MarkerFaceColor', mutColor, 'MarkerSize', 10, 'LineWidth', 2, 'CapSize', 6);

% Add to legend arrays
plotHandles = [h_wt, h_mut];
legendTexts = {'Wild-type (Mean ± SEM)', 'Mutant (Mean ± SEM)'};

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
    
    wtMeans_phase = zeros(1, numDays);
    wtSEMs_phase = zeros(1, numDays);
    mutMeans_phase = zeros(1, numDays);
    mutSEMs_phase = zeros(1, numDays);
    
    % Initialize handles for legend
    plotHandles_phase = [];
    legendTexts_phase = {};
    
    % Plot individual data points for wild-type mice
    for d = 1:numDays
        dayKey = sprintf('day%d', d);
        
        if isfield(dayStats.wild_type, dayKey) && isfield(dayStats.wild_type.(dayKey).(phase), 'counts') && ...
           ~isempty(dayStats.wild_type.(dayKey).(phase).counts)
            individualCounts = dayStats.wild_type.(dayKey).(phase).counts;
            % Add small random jitter to x-position for visibility
            xPos = d + (rand(size(individualCounts)) - 0.5) * 0.1;
            scatter(xPos, individualCounts, 50, wtColor, 'o', 'filled', 'MarkerFaceAlpha', 0.3);
            
            % Store mean and SEM
            wtMeans_phase(d) = dayStats.wild_type.(dayKey).(phase).mean_count;
            wtSEMs_phase(d) = dayStats.wild_type.(dayKey).(phase).sem_count;
        end
    end
    
    % Plot individual data points for mutant mice
    for d = 1:numDays
        dayKey = sprintf('day%d', d);
        
        if isfield(dayStats.mutant, dayKey) && isfield(dayStats.mutant.(dayKey).(phase), 'counts') && ...
           ~isempty(dayStats.mutant.(dayKey).(phase).counts)
            individualCounts = dayStats.mutant.(dayKey).(phase).counts;
            % Add small random jitter to x-position for visibility
            xPos = d + (rand(size(individualCounts)) - 0.5) * 0.1 + 0.2;
            scatter(xPos, individualCounts, 50, mutColor, 'o', 'filled', 'MarkerFaceAlpha', 0.3);
            
            % Store mean and SEM
            mutMeans_phase(d) = dayStats.mutant.(dayKey).(phase).mean_count;
            mutSEMs_phase(d) = dayStats.mutant.(dayKey).(phase).sem_count;
        end
    end
    
    % Plot means with error bars on top
    h_wt_phase = errorbar(dayNumbers, wtMeans_phase, wtSEMs_phase, 'Color', wtColor, 'LineStyle', 'none', ...
        'Marker', 'o', 'MarkerFaceColor', wtColor, 'MarkerSize', 10, 'LineWidth', 2, 'CapSize', 6);
    h_mut_phase = errorbar(dayNumbers + 0.2, mutMeans_phase, mutSEMs_phase, 'Color', mutColor, 'LineStyle', 'none', ...
        'Marker', 'o', 'MarkerFaceColor', mutColor, 'MarkerSize', 10, 'LineWidth', 2, 'CapSize', 6);
    
    % Add to legend arrays
    plotHandles_phase = [h_wt_phase, h_mut_phase];
    legendTexts_phase = {'Wild-type (Mean ± SEM)', 'Mutant (Mean ± SEM)'};
    
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

%% ZT hour analysis
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
ztData.wild_type = struct('counts', zeros(0, numZTHours), 'totalDurations', zeros(0, numZTHours), 'avgDurations', zeros(0, numZTHours));
ztData.mutant = struct('counts', zeros(0, numZTHours), 'totalDurations', zeros(0, numZTHours), 'avgDurations', zeros(0, numZTHours));

% Create structures to store ZT data by day for each mouse
ztDataByDay = struct();
ztDataByDay.wild_type = struct();
ztDataByDay.mutant = struct();

% Process each CSV file for ZT hour analysis
fprintf('Processing files for ZT hour analysis...\n');

for i = 1:length(csvFiles)
    csvPath = fullfile(csvFiles(i).folder, csvFiles(i).name);
    currentMouseID = mouseIDs{i};
    
    % Determine genotype based on user input
    if ismember(currentMouseID, wtMice)
        mouseGenotype = 'wild_type';
    elseif ismember(currentMouseID, mutMice)
        mouseGenotype = 'mutant';
    else
        % Skip mice not assigned to any genotype
        continue;
    end
    
    fprintf('Processing ZT data for mouse %s (genotype: %s)...\n', currentMouseID, strrep(mouseGenotype, '_', '-'));
    
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

for genotype = {'wild_type', 'mutant'}
    gen = genotype{1};
    if isfield(ztData, gen)
        % Calculate statistics for bout counts
        if ~isempty(ztData.(gen).counts)
            ztMean.(gen).counts = mean(ztData.(gen).counts, 1);
            ztStd.(gen).counts = std(ztData.(gen).counts, 0, 1);
            ztSem.(gen).counts = ztStd.(gen).counts / sqrt(size(ztData.(gen).counts, 1));
        else
            fprintf('No ZT bout count data found for genotype: %s\n', gen);
            ztMean.(gen).counts = zeros(1, numZTHours);
            ztStd.(gen).counts = zeros(1, numZTHours);
            ztSem.(gen).counts = zeros(1, numZTHours);
        end
        
        % Calculate statistics for total durations
        if ~isempty(ztData.(gen).totalDurations)
            ztMean.(gen).totalDurations = mean(ztData.(gen).totalDurations, 1);
            ztStd.(gen).totalDurations = std(ztData.(gen).totalDurations, 0, 1);
            ztSem.(gen).totalDurations = ztStd.(gen).totalDurations / sqrt(size(ztData.(gen).totalDurations, 1));
        else
            fprintf('No ZT total duration data found for genotype: %s\n', gen);
            ztMean.(gen).totalDurations = zeros(1, numZTHours);
            ztStd.(gen).totalDurations = zeros(1, numZTHours);
            ztSem.(gen).totalDurations = zeros(1, numZTHours);
        end
        
        % Calculate statistics for average durations
        if ~isempty(ztData.(gen).avgDurations)
            ztMean.(gen).avgDurations = mean(ztData.(gen).avgDurations, 1);
            ztStd.(gen).avgDurations = std(ztData.(gen).avgDurations, 0, 1);
            ztSem.(gen).avgDurations = ztStd.(gen).avgDurations / sqrt(size(ztData.(gen).avgDurations, 1));
        else
            fprintf('No ZT average duration data found for genotype: %s\n', gen);
            ztMean.(gen).avgDurations = zeros(1, numZTHours);
            ztStd.(gen).avgDurations = zeros(1, numZTHours);
            ztSem.(gen).avgDurations = zeros(1, numZTHours);
        end
    end
end

%% Calculate per-day ZT statistics by genotype
fprintf('Calculating per-day ZT statistics...\n');

% Initialize structures for per-day ZT statistics
ztStatsByDay = struct();
ztStatsByDay.wild_type = struct();
ztStatsByDay.mutant = struct();

for d = 1:numDays
    dayKey = sprintf('day%d', d);
    
    for genotype = {'wild_type', 'mutant'}
        gen = genotype{1};
        
        % Initialize storage for this day and genotype
        ztStatsByDay.(gen).(dayKey) = struct();
        ztStatsByDay.(gen).(dayKey).counts = struct('data', [], 'mean', zeros(1, numZTHours), 'std', zeros(1, numZTHours), 'sem', zeros(1, numZTHours));
        ztStatsByDay.(gen).(dayKey).totalDurations = struct('data', [], 'mean', zeros(1, numZTHours), 'std', zeros(1, numZTHours), 'sem', zeros(1, numZTHours));
        ztStatsByDay.(gen).(dayKey).avgDurations = struct('data', [], 'mean', zeros(1, numZTHours), 'std', zeros(1, numZTHours), 'sem', zeros(1, numZTHours));
        
        % Collect data from all mice of this genotype for this day
        if isfield(ztDataByDay, gen)
            mouseFields = fieldnames(ztDataByDay.(gen));
            
            for m = 1:length(mouseFields)
                mouseFieldName = mouseFields{m};
                mouseData = ztDataByDay.(gen).(mouseFieldName);
                
                % Add this mouse's data for this day
                ztStatsByDay.(gen).(dayKey).counts.data = [ztStatsByDay.(gen).(dayKey).counts.data; mouseData.counts(d, :)];
                ztStatsByDay.(gen).(dayKey).totalDurations.data = [ztStatsByDay.(gen).(dayKey).totalDurations.data; mouseData.totalDurations(d, :)];
                ztStatsByDay.(gen).(dayKey).avgDurations.data = [ztStatsByDay.(gen).(dayKey).avgDurations.data; mouseData.avgDurations(d, :)];
            end
            
            % Calculate statistics if we have data
            if ~isempty(ztStatsByDay.(gen).(dayKey).counts.data)
                ztStatsByDay.(gen).(dayKey).counts.mean = mean(ztStatsByDay.(gen).(dayKey).counts.data, 1);
                ztStatsByDay.(gen).(dayKey).counts.std = std(ztStatsByDay.(gen).(dayKey).counts.data, 0, 1);
                ztStatsByDay.(gen).(dayKey).counts.sem = ztStatsByDay.(gen).(dayKey).counts.std / sqrt(size(ztStatsByDay.(gen).(dayKey).counts.data, 1));
                
                ztStatsByDay.(gen).(dayKey).totalDurations.mean = mean(ztStatsByDay.(gen).(dayKey).totalDurations.data, 1);
                ztStatsByDay.(gen).(dayKey).totalDurations.std = std(ztStatsByDay.(gen).(dayKey).totalDurations.data, 0, 1);
                ztStatsByDay.(gen).(dayKey).totalDurations.sem = ztStatsByDay.(gen).(dayKey).totalDurations.std / sqrt(size(ztStatsByDay.(gen).(dayKey).totalDurations.data, 1));
                
                ztStatsByDay.(gen).(dayKey).avgDurations.mean = mean(ztStatsByDay.(gen).(dayKey).avgDurations.data, 1);
                ztStatsByDay.(gen).(dayKey).avgDurations.std = std(ztStatsByDay.(gen).(dayKey).avgDurations.data, 0, 1);
                ztStatsByDay.(gen).(dayKey).avgDurations.sem = ztStatsByDay.(gen).(dayKey).avgDurations.std / sqrt(size(ztStatsByDay.(gen).(dayKey).avgDurations.data, 1));
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
for genotype = {'wild_type', 'mutant'}
    gen = genotype{1};
    
    if isfield(animalIDs, gen) && ~isempty(animalIDs.(gen))
        for m = 1:min(length(animalIDs.(gen)), size(ztData.(gen).counts, 1))
            mouseID = animalIDs.(gen){m};
            
            % Add counts data
            ztCountsTable.([mouseID '_Counts']) = ztData.(gen).counts(m,:)';
            
            % Add total durations data
            ztTotalDurationsTable.([mouseID '_Total']) = ztData.(gen).totalDurations(m,:)';
            
            % Add average durations data
            ztAvgDurationsTable.([mouseID '_Avg']) = ztData.(gen).avgDurations(m,:)';
        end
    end
end

% Add genotype statistics to tables
for genotype = {'wild_type', 'mutant'}
    gen = genotype{1};
    display_name = strrep(gen, '_', '-'); % Convert 'wild_type' to 'wild-type'
    
    if isfield(ztMean, gen)
        % Add bout counts statistics
        ztCountsTable.([display_name '_Mean']) = ztMean.(gen).counts';
        ztCountsTable.([display_name '_SD']) = ztStd.(gen).counts';
        ztCountsTable.([display_name '_SEM']) = ztSem.(gen).counts';
        
        % Add total durations statistics
        ztTotalDurationsTable.([display_name '_Mean']) = ztMean.(gen).totalDurations';
        ztTotalDurationsTable.([display_name '_SD']) = ztStd.(gen).totalDurations';
        ztTotalDurationsTable.([display_name '_SEM']) = ztSem.(gen).totalDurations';
        
        % Add average durations statistics
        ztAvgDurationsTable.([display_name '_Mean']) = ztMean.(gen).avgDurations';
        ztAvgDurationsTable.([display_name '_SD']) = ztStd.(gen).avgDurations';
        ztAvgDurationsTable.([display_name '_SEM']) = ztSem.(gen).avgDurations';
    end
end

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
    if isfield(animalIDs, 'wild_type')
        allMouseIDs_ZT = [allMouseIDs_ZT, animalIDs.wild_type];
    end
    if isfield(animalIDs, 'mutant')
        allMouseIDs_ZT = [allMouseIDs_ZT, animalIDs.mutant];
    end
    
    for m = 1:length(allMouseIDs_ZT)
        mouseID = allMouseIDs_ZT{m};
        mouseFieldName = makeValidFieldName(mouseID);
        
        % Find which genotype this mouse belongs to
        if isfield(animalIDs, 'wild_type') && any(strcmp(animalIDs.wild_type, mouseID))
            gen = 'wild_type';
        elseif isfield(animalIDs, 'mutant') && any(strcmp(animalIDs.mutant, mouseID))
            gen = 'mutant';
        else
            continue;
        end
        
        % Check if this mouse has ZT data
        if isfield(ztDataByDay, gen) && isfield(ztDataByDay.(gen), mouseFieldName)
            mouseZTData = ztDataByDay.(gen).(mouseFieldName);
            
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
    
    % Add genotype statistics for this day
    for genotype = {'wild_type', 'mutant'}
        gen = genotype{1};
        display_name = strrep(gen, '_', '-');
        
        if isfield(ztStatsByDay, gen) && isfield(ztStatsByDay.(gen), dayKey)
            dayStats = ztStatsByDay.(gen).(dayKey);
            
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

%% Create ZT hour plots
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
    
    % Plot individual data points for wild-type
    if isfield(ztData, 'wild_type') && ~isempty(ztData.wild_type.(plotType))
        for i = 1:size(ztData.wild_type.(plotType), 1)
            scatter(0:23, ztData.wild_type.(plotType)(i,:), 50, wtColor, 'o', 'filled', 'MarkerFaceAlpha', 0.3);
        end
        
        % Plot mean with error bars and capture handle
        h_wt = errorbar(0:23, ztMean.wild_type.(plotType), ztSem.wild_type.(plotType), ...
            'Color', wtColor, 'LineWidth', 2, 'Marker', 'o', 'MarkerFaceColor', wtColor, 'MarkerSize', 10);
        
        % Add to legend arrays
        plotHandles = [plotHandles, h_wt];
        legendTexts{end+1} = 'Wild-type (Mean ± SEM)';
    end
    
    % Plot individual data points for mutant
    if isfield(ztData, 'mutant') && ~isempty(ztData.mutant.(plotType))
        for i = 1:size(ztData.mutant.(plotType), 1)
            scatter((0:23)+0.2, ztData.mutant.(plotType)(i,:), 50, mutColor, 'o', 'filled', 'MarkerFaceAlpha', 0.3);
        end
        
        % Plot mean with error bars and capture handle
        h_mut = errorbar((0:23)+0.2, ztMean.mutant.(plotType), ztSem.mutant.(plotType), ...
            'Color', mutColor, 'LineWidth', 2, 'Marker', 'o', 'MarkerFaceColor', mutColor, 'MarkerSize', 10);
        
        % Add to legend arrays
        plotHandles = [plotHandles, h_mut];
        legendTexts{end+1} = 'Mutant (Mean ± SEM)';
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
dayColors = colormap(lines(numDays)); % Generates distinct colors for each day

% 1. ZT bout counts across days by genotype
for genotype = {'wild_type', 'mutant'}
    gen = genotype{1};
    display_name = strrep(gen, '_', '-');
    
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
        
        if isfield(ztStatsByDay, gen) && isfield(ztStatsByDay.(gen), dayKey) && ...
           ~isempty(ztStatsByDay.(gen).(dayKey).counts.mean)
            
            meanData = ztStatsByDay.(gen).(dayKey).counts.mean;
            semData = ztStatsByDay.(gen).(dayKey).counts.sem;
            
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
    saveas(gcf, fullfile(outputFolder, ['ZT_Counts_Across_Days_' display_name '.fig']));
    saveas(gcf, fullfile(outputFolder, ['ZT_Counts_Across_Days_' display_name '.png']));
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
    
    % Plot both genotypes
    plotHandles = [];
    legendTexts = {};
    
    for genotype = {'wild_type', 'mutant'}
        gen = genotype{1};
        display_name = strrep(gen, '_', '-');
        
        if strcmp(gen, 'wild_type')
            color = wtColor;
        else
            color = mutColor;
        end
        
        if isfield(ztStatsByDay, gen) && isfield(ztStatsByDay.(gen), dayKey) && ...
           ~isempty(ztStatsByDay.(gen).(dayKey).counts.mean)
            
            meanData = ztStatsByDay.(gen).(dayKey).counts.mean;
            semData = ztStatsByDay.(gen).(dayKey).counts.sem;
            
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
if isfield(animalIDs, 'wild_type')
    allMouseIDs_plot = [allMouseIDs_plot, animalIDs.wild_type];
end
if isfield(animalIDs, 'mutant')
    allMouseIDs_plot = [allMouseIDs_plot, animalIDs.mutant];
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
        if isfield(animalIDs, 'wild_type') && any(strcmp(animalIDs.wild_type, mouseID))
            gen = 'wild_type';
            mouseColor = wtColor;
        elseif isfield(animalIDs, 'mutant') && any(strcmp(animalIDs.mutant, mouseID))
            gen = 'mutant';
            mouseColor = mutColor;
        else
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
        if isfield(ztDataByDay, gen) && isfield(ztDataByDay.(gen), mouseFieldName)
            mouseZTData = ztDataByDay.(gen).(mouseFieldName);
            
            for d = 1:numDays
                plot(0:23, mouseZTData.counts(d, :), ...
                    'Color', dayColors(d, :), 'LineWidth', 1.5, 'Marker', 'o', ...
                    'MarkerFaceColor', dayColors(d, :), 'MarkerSize', 4);
            end
        end
        
        % Bring plot lines to front
        uistack(findobj(gca, 'Type', 'line'), 'top');
        
        title(['Mouse ' mouseID ' (' strrep(gen, '_', '-') ')']);
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

%% Create ZT plots across experimental days
fprintf('Creating ZT bout counts across experimental days plots...\n');

% 1. ZT bout counts across days by genotype
for genotype = {'wild_type', 'mutant'}
    gen = genotype{1};
    display_name = strrep(gen, '_', '-');
    
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
        
        if isfield(ztStatsByDay, gen) && isfield(ztStatsByDay.(gen), dayKey) && ...
           ~isempty(ztStatsByDay.(gen).(dayKey).counts.mean)
            
            meanData = ztStatsByDay.(gen).(dayKey).counts.mean;
            semData = ztStatsByDay.(gen).(dayKey).counts.sem;
            
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
    saveas(gcf, fullfile(outputFolder, ['ZT_Counts_Across_Days_' display_name '.fig']));
    saveas(gcf, fullfile(outputFolder, ['ZT_Counts_Across_Days_' display_name '.png']));
end

% Restore normal figure visibility
set(0, 'DefaultFigureVisible', 'on');

fprintf('ZT across days analysis complete!\n');
fprintf('All figures saved to %s\n', outputFolder);
fprintf('Analysis complete!\n');