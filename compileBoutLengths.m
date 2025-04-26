% Script to analyze mouse sleep data by genotype from CSV files
% This script:
% 1. Finds all CSV files in input folder (no subfolders)
% 2. Groups mice by genotype (wild-type or mutant)
% 3. Retrieves epoch durations for each mouse from CSV files
% 4. Sorts sleep bouts into light phase (6:00-18:00) and dark phase (18:00-6:00)
% 5. Bins the epochs by their duration (0s, 4s, 8s, 16s, 32s, 64s, 128s, 256s, 512s)
% 6. Computes averages and standard deviations by genotype (overall, light phase, dark phase)
% 7. Plots results as line, bar, and dot graphs with error bars for all data categories

%% Parameters and Setup
mainFolder = '/Users/davidrivas/Documents/research/sleep-box/09-18-23'; % Use current directory, or specify your path
genotypes = {'wild-type', 'mutant'}; % Define genotypes 

% Define bin edges for epoch durations (in seconds)
binEdges = [0, 4, 8, 16, 32, 64, 128, 256, 512, inf];
binLabels = {'0-4s', '4-8s', '8-16s', '16-32s', '32-64s', '64-128s', '128-256s', '256-512s', '>512s'};
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
                            
                            % Determine phase based on hour (6:00-18:00 is light phase)
                            if hour >= 6 && hour < 18
                                % Light phase
                                lightDurations = [lightDurations; sleepDurations(t)];
                            else
                                % Dark phase
                                darkDurations = [darkDurations; sleepDurations(t)];
                            end
                        else
                            fprintf('    Warning: Unable to parse timestamp: %s\n', timestamp);
                        end
                    catch e
                        fprintf('    Warning: Error processing timestamp %s: %s\n', timestamp, e.message);
                    end
                end
                
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

% 1. Create table for individual mouse data (all sleep)
mouseDataTable = array2table(zeros(length(allBinLabels), 0));
mouseDataTable.Bin = allBinLabels;

fprintf('Preparing individual mouse data sheet (all sleep)...\n');
% Add columns for each individual mouse - all sleep data
for genotype = {'wild_type', 'mutant'}
    gen = genotype{1};
    
    if isfield(animalIDs, gen) && ~isempty(animalIDs.(gen))
        for m = 1:length(animalIDs.(gen))
            mouseID = animalIDs.(gen){m};
            fprintf('  Adding all sleep data for mouse: %s\n', mouseID);
            
            % Create data column for all sleep
            if isfield(totalCounts, gen) && length(totalCounts.(gen).sleep) >= m
                sleepData = [totalCounts.(gen).sleep(m); epochCounts.(gen).sleep(m,:)'];
                
                % Add as full column
                mouseDataTable.([mouseID '_Sleep']) = sleepData;
            end
        end
    end
end

% 2. Create table for individual mouse data (light/dark phases)
mouseDataLD_Table = array2table(zeros(length(allBinLabels), 0));
mouseDataLD_Table.Bin = allBinLabels;

fprintf('Preparing individual mouse data sheet (light/dark phases)...\n');
% Add columns for each individual mouse - light/dark phase data
for genotype = {'wild_type', 'mutant'}
    gen = genotype{1};
    
    if isfield(animalIDs, gen) && ~isempty(animalIDs.(gen))
        for m = 1:length(animalIDs.(gen))
            mouseID = animalIDs.(gen){m};
            fprintf('  Adding light/dark phase data for mouse: %s\n', mouseID);
            
            % Create data columns for light and dark phases
            if isfield(totalCounts, gen) && length(totalCounts.(gen).light) >= m
                lightData = [totalCounts.(gen).light(m); epochCounts.(gen).light(m,:)'];
                darkData = [totalCounts.(gen).dark(m); epochCounts.(gen).dark(m,:)'];
                
                % Add as full columns
                mouseDataLD_Table.([mouseID '_Light']) = lightData;
                mouseDataLD_Table.([mouseID '_Dark']) = darkData;
            end
        end
    end
end

% 3. Create table for genotype statistics (all sleep)
genotypeStatsTable = array2table(zeros(length(allBinLabels), 0));
genotypeStatsTable.Bin = allBinLabels;

fprintf('Preparing genotype statistics sheet (all sleep)...\n');
% Add columns for each genotype with simplified names - all sleep data
for genotype = {'wild_type', 'mutant'}
    gen = genotype{1};
    display_name = strrep(gen, '_', '-'); % Convert 'wild_type' to 'wild-type'
    
    if isfield(totalCounts, gen)
        % Mean columns - All Sleep
        genotypeStatsTable.([display_name '_Sleep_Mean']) = [mean(totalCounts.(gen).sleep); meanCounts.(gen).sleep'];
        
        % SD columns - All Sleep
        genotypeStatsTable.([display_name '_Sleep_SD']) = [std(totalCounts.(gen).sleep); stdCounts.(gen).sleep'];
        
        % SEM columns - All Sleep
        totalSEM_sleep = std(totalCounts.(gen).sleep)/sqrt(length(totalCounts.(gen).sleep));
        genotypeStatsTable.([display_name '_Sleep_SEM']) = [totalSEM_sleep; semCounts.(gen).sleep'];
    end
end

% 4. Create table for genotype statistics (light/dark phases)
genotypeStatsLD_Table = array2table(zeros(length(allBinLabels), 0));
genotypeStatsLD_Table.Bin = allBinLabels;

fprintf('Preparing genotype statistics sheet (light/dark phases)...\n');
% Add columns for each genotype with simplified names - light/dark phase data
for genotype = {'wild_type', 'mutant'}
    gen = genotype{1};
    display_name = strrep(gen, '_', '-'); % Convert 'wild_type' to 'wild-type'
    
    if isfield(totalCounts, gen)
        % Mean columns - Light Phase
        genotypeStatsLD_Table.([display_name '_Light_Mean']) = [mean(totalCounts.(gen).light); meanCounts.(gen).light'];
        
        % SD columns - Light Phase
        genotypeStatsLD_Table.([display_name '_Light_SD']) = [std(totalCounts.(gen).light); stdCounts.(gen).light'];
        
        % SEM columns - Light Phase
        totalSEM_light = std(totalCounts.(gen).light)/sqrt(length(totalCounts.(gen).light));
        genotypeStatsLD_Table.([display_name '_Light_SEM']) = [totalSEM_light; semCounts.(gen).light'];
        
        % Mean columns - Dark Phase
        genotypeStatsLD_Table.([display_name '_Dark_Mean']) = [mean(totalCounts.(gen).dark); meanCounts.(gen).dark'];
        
        % SD columns - Dark Phase
        genotypeStatsLD_Table.([display_name '_Dark_SD']) = [std(totalCounts.(gen).dark); stdCounts.(gen).dark'];
        
        % SEM columns - Dark Phase
        totalSEM_dark = std(totalCounts.(gen).dark)/sqrt(length(totalCounts.(gen).dark));
        genotypeStatsLD_Table.([display_name '_Dark_SEM']) = [totalSEM_dark; semCounts.(gen).dark'];
    end
end

% 5. Create a summary table mapping mice to genotypes
mouseGenotypeSummary = table();
mouseGenotypeSummary.MouseID = [animalIDs.wild_type(:); animalIDs.mutant(:)];
genotypes_list = [repmat({'wild-type'}, length(animalIDs.wild_type), 1); 
                  repmat({'mutant'}, length(animalIDs.mutant), 1)];
mouseGenotypeSummary.Genotype = genotypes_list;

% Write all tables to different sheets in the Excel file
xlsFilePath = fullfile(mainFolder, 'sleep_bout_analysis.xls');
fprintf('Writing data to Excel file: %s\n', xlsFilePath);

% Write to Excel file - all five sheets
writetable(mouseDataTable, xlsFilePath, 'Sheet', 'Mouse_Data');
writetable(mouseDataLD_Table, xlsFilePath, 'Sheet', 'Mouse_Data_LD');
writetable(genotypeStatsTable, xlsFilePath, 'Sheet', 'Genotype_Statistics');
writetable(genotypeStatsLD_Table, xlsFilePath, 'Sheet', 'Genotype_Statistics_LD');
writetable(mouseGenotypeSummary, xlsFilePath, 'Sheet', 'Mouse_Genotypes');

fprintf('Excel file created successfully.\n');

%% Create figures
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
    
    %% 1. Line Plot
    figure('Name', [categoryLabel ' Bout Durations - Line Plot'], 'Position', [100, 100, 800, 500]);
    hold on;
    
    % Plot wild-type data with error bars
    if isfield(meanCounts, 'wild_type') && all(isfinite(meanCounts.wild_type.(category)))
        errorbar(1:numBins, meanCounts.wild_type.(category), semCounts.wild_type.(category), ...
            'Color', wtColor, 'LineWidth', 2, 'Marker', 'o', 'MarkerFaceColor', wtColor, 'MarkerSize', 8);
    end
    
    % Plot mutant data with error bars
    if isfield(meanCounts, 'mutant') && all(isfinite(meanCounts.mutant.(category)))
        errorbar(1:numBins, meanCounts.mutant.(category), semCounts.mutant.(category), ...
            'Color', mutColor, 'LineWidth', 2, 'Marker', 'o', 'MarkerFaceColor', mutColor, 'MarkerSize', 8);
    end
    
    % Add labels and legend
    title([categoryLabel ' Bout Durations']);
    xlabel('Bout Duration');
    ylabel('Number of Bouts');
    set(gca, 'XTick', 1:numBins, 'XTickLabel', binLabels);
    xtickangle(45);
    legend('Wild-type (Mean ± SEM)', 'Mutant (Mean ± SEM)', 'Location', 'best');
    grid on;
    hold off;
    
    % Save the figure
    saveas(gcf, fullfile(outputFolder, [strrep(categoryLabel, ' ', '_') '_Bouts_Line.fig']));
    saveas(gcf, fullfile(outputFolder, [strrep(categoryLabel, ' ', '_') '_Bouts_Line.png']));
    
    %% 2. Bar Plot
    figure('Name', [categoryLabel ' Bout Durations - Bar Plot'], 'Position', [100, 100, 800, 500]);
    hold on;
    
    % Bar width and positions
    barWidth = 0.35;
    wtPos = (1:numBins) - barWidth/2;
    mutPos = (1:numBins) + barWidth/2;
    
    % Clear barHandles for this new figure
    clear barHandles;
    
    % Plot wild-type bars
    if isfield(meanCounts, 'wild_type') && all(isfinite(meanCounts.wild_type.(category)))
        barHandles(1) = bar(wtPos, meanCounts.wild_type.(category), barWidth, 'FaceColor', wtColor);
        
        % Add error bars
        errorbar(wtPos, meanCounts.wild_type.(category), stdCounts.wild_type.(category), '.k');
    end
    
    % Plot mutant bars
    if isfield(meanCounts, 'mutant') && all(isfinite(meanCounts.mutant.(category)))
        barHandles(2) = bar(mutPos, meanCounts.mutant.(category), barWidth, 'FaceColor', mutColor);
        
        % Add error bars
        errorbar(mutPos, meanCounts.mutant.(category), stdCounts.mutant.(category), '.k');
    end
    
    % Add labels and legend
    title([categoryLabel ' Bout Durations']);
    xlabel('Bout Duration');
    ylabel('Number of Bouts');
    set(gca, 'XTick', 1:numBins, 'XTickLabel', binLabels);
    xtickangle(45);
    
    % Add legend if both genotypes have data
    if exist('barHandles', 'var') && length(barHandles) >= 2
        legend(barHandles, {'Wild-type (Mean ± SD)', 'Mutant (Mean ± SD)'}, 'Location', 'best');
    end
    
    grid on;
    hold off;
    
    % Save the figure
    saveas(gcf, fullfile(outputFolder, [strrep(categoryLabel, ' ', '_') '_Bouts_Bar.fig']));
    saveas(gcf, fullfile(outputFolder, [strrep(categoryLabel, ' ', '_') '_Bouts_Bar.png']));
    
    %% 3. Dot Plot (Individual mice with mean)
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
fprintf('Analysis complete!\n');