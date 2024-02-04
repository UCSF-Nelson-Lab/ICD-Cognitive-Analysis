%%
clear all;
close all;
% Specify the Excel file path
filePath = '{Your file path}';
% Get the sheet names in the Excel file
[~, sheetNames] = xlsfinfo(excelFilePath);
% Initialize cell arrays to store averages and standard errors for each sheet
averages = cell(numel(sheetNames), 1);
stdErrors = cell(numel(sheetNames), 1);
% Define the column index for "Animal ID" (assuming it's in the first column)
animalIDColumn = 1;

%% calculate the averages across different days of PPX injection
% Initialize cell arrays to store averages and standard errors for each sheet
SeparateAvg = cell(numel(sheetNames), 1);
SeparateSTE = cell(numel(sheetNames), 1);

% Define the column index for "Animal ID" (assuming it's in the first column)
animalIDColumn = 1;
% Define the list of animal IDs to exclude
excludedAnimalIDs = {''};
% Loop through each sheet
for sheetIndex = 1:numel(sheetNames)
    % Read the data for the current sheet
    [~, ~, raw] = xlsread(excelFilePath, sheetNames{sheetIndex});
    % Initialize a structure to hold the rows for each unique Animal ID
    rowsData = struct('row1', [], 'row2', [], 'row3', [], 'row4', []);
    
    % Get the unique animal IDs
    rawAnimalIDs = raw(2:end, animalIDColumn);  % Assuming 1st row is header
    uniqueAnimalIDs = unique(rawAnimalIDs);
    
    % Loop through each unique Animal ID
    for i = 1:length(uniqueAnimalIDs)
        if ismember(uniqueAnimalIDs{i}, excludedAnimalIDs)
            continue;
        end
        id = uniqueAnimalIDs{i};        
        % Find the rows that have this unique Animal ID
        rowIdx = find(strcmp(rawAnimalIDs, id));
       
        % Extract up to first 4 rows for each Animal ID
        for j = 1:min(4, length(rowIdx))
            if rowIdx(j) + 1 <= size(raw, 1)  % +1 because 1st row is header
                rowNumeric = cellfun(@(x) iff(isnumeric(x), x, NaN), raw(rowIdx(j) + 1, 2:end), 'UniformOutput', false);
                rowsData.(sprintf('row%d', j)) = [rowsData.(sprintf('row%d', j)); cell2mat(rowNumeric)];
            end
        end
    end   
    % Calculate averages and standard errors for each set of rows
    for j = 1:4
        fieldName = sprintf('row%d', j);
        colAverages = nanmean(rowsData.(fieldName), 1);
        colStdErrors = nanstd(rowsData.(fieldName), 0, 1) ./ sqrt(size(rowsData.(fieldName), 1));       
        % Store the results in cell arrays for this sheet
        SeparateAvg{sheetIndex}.(fieldName) = colAverages;
        SeparateSTE{sheetIndex}.(fieldName) = colStdErrors;
     end
end

%% plot the delay discounting curve across days of PPX injections
% Define the custom titles for the plots
customTitles = {'park-baseline', 'park-PPX', 'park-PPXwashout', 'saline-baseline', 'saline-PPX', 'saline-PPX washout'};
% Initialize empty arrays to store plotting data for all structs
allDataToPlot = [];
allStdErrToPlot = [];
allStructNames = {};
allRowNames = {};
delay = [0, 3, 6, 9];
columnsToPlot = [6, 15, 24, 35];
% % # of LR choice:[6, 15, 24, 35]
  % # of center omission:[7, 16, 25, 36]
  % # of side omission: [8, 17, 26, 37]
  % reaction time: [9, 18, 27, 38]
  % trigger time for LR: [10, 19, 28, 39]
  % trigger time for SR: [11, 20, 29, 40]
  % reaction trigger difference for L/SR: [12, 21, 30, 41]

% Loop through each struct (from struct1 to struct6)
% for structIndex = 1:6
for structIndex = 7:12
    % Retrieve the current struct from SeparateAvg and SeparateSTE cell arrays
    currAvgStruct = SeparateAvg{structIndex};
    currStdStruct = SeparateSTE{structIndex};
    
    % Initialize arrays for the data to plot
    dataToPlot = zeros(4, numel(columnsToPlot));
    stdErrToPlot = zeros(4, numel(columnsToPlot));
    
    % Extract data for the specified columns from each row
    for i = 1:4
        rowName = sprintf('row%d', i);
        dataToPlot(i, :) = currAvgStruct.(rowName)(columnsToPlot);
        stdErrToPlot(i, :) = currStdStruct.(rowName)(columnsToPlot);
    end
    
    % Store plotting data for current struct in the global arrays
    allDataToPlot = [allDataToPlot; dataToPlot];
    allStdErrToPlot = [allStdErrToPlot; stdErrToPlot];
    
    % Create row names for the current struct
    currentRowNames = arrayfun(@(i) sprintf('Struct%d_Row%d', structIndex, i), 1:4, 'UniformOutput', false);
    allRowNames = [allRowNames, currentRowNames];
    
    % Store the struct name for labeling columns later
    allStructNames = [allStructNames, customTitles(structIndex)];
    
    % Create a new figure for each struct
    figure;
    % Plotting with error bars
    hold on;
    for i = 7:9
        errorbar(delay, dataToPlot(i, :), stdErrToPlot(i, :), '-o', 'DisplayName', sprintf('Row %d', i));
    end
    hold off;
    
    % Add labels and legend
    xlabel('Delay');
    ylabel('Average Value');
    title(customTitles{structIndex});
    legend('show');
end
% Create tables for easy visualization
allDataToPlotTable = array2table(allDataToPlot, 'RowNames', allRowNames, 'VariableNames', arrayfun(@(x) ['Delay_' num2str(x)], delay, 'UniformOutput', false));
allStdErrToPlotTable = array2table(allStdErrToPlot, 'RowNames', allRowNames, 'VariableNames', arrayfun(@(x) ['Delay_' num2str(x)], delay, 'UniformOutput', false));

%% calculate the averages and standard error for each individual mice
sheetAverages = cell(numel(sheetNames), 1); % Initialize to store averages
extractedColumnsForSheet = cell(numel(sheetNames), 1); % Initialize to store selected columns

% Loop through each sheet for averaging
for sheetIndex = 1:numel(sheetNames)
    [~, ~, raw] = xlsread(excelFilePath, sheetNames{sheetIndex});
    animalIDs = raw(2:end, 1);  % Assuming animal IDs are in the first column and the first row is a header
    uniqueAnimalIDs = unique(animalIDs);   
    % Initialization
    numUniqueAnimals = numel(uniqueAnimalIDs);
    averagesForSheet = nan(numUniqueAnimals, size(raw, 2) - 1);

    for i = 1:numUniqueAnimals
        id = uniqueAnimalIDs{i};
        rows = find(strcmp(animalIDs, id));
        rows = rows + 1;  % Adjusting index for 1-based indexing
        animalData = raw(rows, 2:end);  % Assuming data starts from the 2nd column
        
        % Convert empty cells to NaN
        animalData(cellfun(@isempty, animalData)) = {NaN};
        animalData(cellfun(@(x) ~isnumeric(x) && ~islogical(x), animalData)) = {NaN};
        
        % Convert to numeric matrix
        animalDataNumeric = cell2mat(animalData);
        
        % Column-wise averages
        averagesForSheet(i, :) = nanmean(animalDataNumeric, 1);
    end
    
    % Store the averages for this sheet
    sheetAverages{sheetIndex} = averagesForSheet;
    
    % Extract columns 6, 15, 24, 35
    extractedColumns = averagesForSheet(:, [6, 15, 24, 35]);
    
    % Store the extracted columns for this sheet
    extractedColumnsForSheet{sheetIndex} = extractedColumns;
end


%% plot the overall delay discounting curve/other parameters from each condition
delay = [0, 3, 6, 9];
cols_to_average = [6, 15, 24, 35];
% % # of LR choice:[6, 15, 24, 35]
  % # of center omission:[7, 16, 25, 36]
  % # of side omission: [8, 17, 26, 37]
  % reaction time: [9, 18, 27, 38]
  % trigger time for LR: [10, 19, 28, 39]
  % trigger time for SR: [11, 20, 29, 40]
  % reaction trigger difference for L/SR: [12, 21, 30, 41]
  
% Initialize matrices to hold the column-wise averages and standard errors for each sheet
sheetColAverages = cell(numel(sheetNames), 1);
sheetColStdErrors = cell(numel(sheetNames), 1);
% Loop through each sheet
for sheetIndex = 1:numel(sheetNames)
    % Extract the data for the current sheet
    dataForSheet = sheetAverages{sheetIndex};   
    % Calculate the column-wise averages and standard errors for the specified columns
    averages = nanmean(dataForSheet(:, cols_to_average), 1);
    std_errors = nanstd(dataForSheet(:, cols_to_average), 0, 1) ./ sqrt(sum(~isnan(dataForSheet(:, cols_to_average))));   
    % Store the column-wise averages and standard errors for the current sheet
    sheetColAverages{sheetIndex} = averages;
    sheetColStdErrors{sheetIndex} = std_errors;  
end
% Now, plot these averages and standard errors as a function of delay
figure;
% First subplot for the first three cells
subplot(2, 1, 1);
hold on;
for sheetIndex = 7:8  % Assume the first 3 cells are in the first 3 sheets
    % Extract the column-wise averages and standard errors for the current sheet
    averages = sheetColAverages{sheetIndex};
    std_errors = sheetColStdErrors{sheetIndex};
    % Plot the averages as a function of delay
    errorbar(delay, averages, std_errors, '-o', 'DisplayName', ['Sheet ' num2str(sheetIndex)]);
    
end

% Customize the first subplot
xlabel('Delay');
ylabel('%LR choice');
legend('show');
title('Parkinsonian mice');
hold off;

% Second subplot for the last three cells
subplot(2, 1, 2);
hold on;
for sheetIndex = 10:11 % Assume the last 3 cells are in the last sheets numel(sheetNames)
    % Extract the column-wise averages and standard errors for the current sheet
    averages = sheetColAverages{sheetIndex};
    std_errors = sheetColStdErrors{sheetIndex};
    % Plot the averages as a function of delay
    errorbar(delay, averages, std_errors, '-o', 'DisplayName', ['Sheet ' num2str(sheetIndex)]);
    
end

% Customize the second subplot
xlabel('Delay');
ylabel('%LR choice');
legend('show');
title('Healthy mice');
hold off;

%% estimated K, A and AUC values from the fitted delay discounting curve
% Given data
delay = [0, 3, 6, 9];
columnsToFit = [6, 15, 24, 35];

estimatedAValues = cell(numel(sheetNames), 1);
estimatedKValues = cell(numel(sheetNames), 1);
estimatedAUC = cell(numel(sheetNames), 1);
xValuesAtYHalf = cell(numel(sheetNames), 1); % For x values where y = 0.5

for sheetIndex = 1:numel(sheetNames)
    KAdata = sheetAverages{sheetIndex}(:, columnsToFit);    
    AValues = nan(size(KAdata, 1), 1);
    KValues = nan(size(KAdata, 1), 1);
    AUCValues = nan(size(KAdata, 1), 1);
    xAtYHalfValues = nan(size(KAdata, 1), 1); % Initialize x values at y = 0.5
    
    for rowIndex = 1:size(KAdata, 1)
        if all(~isnan(KAdata(rowIndex, :)))
            modelFunc = @(params, x) params(1) ./ (1 + params(2) * x);
            initialParams = [1, 0.5];
            [params, ~, ~, exitflag] = lsqcurvefit(modelFunc, initialParams, delay, KAdata(rowIndex, :), [], [], optimoptions('lsqcurvefit', 'MaxFunEvals', 1000, 'MaxIter', 1000));           
            if exitflag > 0
                AValues(rowIndex) = params(1);
                KValues(rowIndex) = params(2);
                % Calculate AUC
                predictedY = modelFunc(params, delay);
                AUCValues(rowIndex) = trapz(delay, predictedY);
                % Find x where y = 0.5
                xAtYHalf = findXForYHalf(params);
                if xAtYHalf >= 0 && xAtYHalf <= 9
                    xAtYHalfValues(rowIndex) = xAtYHalf;
                else
                    xAtYHalfValues(rowIndex) = NaN;
                end
            else
                AValues(rowIndex) = NaN;
                KValues(rowIndex) = NaN;
                AUCValues(rowIndex) = NaN;
            end
        else
            AValues(rowIndex) = NaN;
            KValues(rowIndex) = NaN;
            AUCValues(rowIndex) = NaN;
        end
    end
    
    estimatedAValues{sheetIndex} = AValues;
    estimatedKValues{sheetIndex} = KValues;
    estimatedAUC{sheetIndex} = AUCValues;
    xValuesAtYHalf{sheetIndex} = xAtYHalfValues; % Store x values for y = 0.5
end

% Plotting
% Estimated A values
figure;
hold on;
firstThreeSheetsA = estimatedAValues(1:3);
minRowsFirstThreeA = min(cellfun(@(sheet) size(sheet, 1), firstThreeSheetsA));
for rowIdx = 1:minRowsFirstThreeA
    plotDataA = cellfun(@(sheet) sheet(rowIdx), firstThreeSheetsA);
    plot(1:3, plotDataA, '-o', 'DisplayName', sprintf('Row %d (A)', rowIdx));
end
xlabel('Sheet Index');
ylabel('Estimated A Value');
title('Estimated A values for First Three Sheets');

% Estimated K values
figure;
hold on;
firstThreeSheetsK = estimatedKValues(1:3);
minRowsFirstThreeK = min(cellfun(@(sheet) size(sheet, 1), firstThreeSheetsK));
for rowIdx = 1:minRowsFirstThreeK
    plotDataK = cellfun(@(sheet) sheet(rowIdx), firstThreeSheetsK);
    plot(1:3, plotDataK, '-o', 'DisplayName', sprintf('Row %d (K)', rowIdx));
end
xlabel('Sheet Index');
ylabel('Estimated K Value');
title('Estimated K values for First Three Sheets');

% Estimated AUC values
figure;
hold on;
firstThreeSheetsAUC = estimatedAUC(1:3);
minRowsFirstThreeAUC = min(cellfun(@(sheet) size(sheet, 1), firstThreeSheetsAUC));
for rowIdx = 1:minRowsFirstThreeAUC
    plotDataAUC = cellfun(@(sheet) sheet(rowIdx), firstThreeSheetsAUC);
    plot(1:3, plotDataAUC, '-o', 'DisplayName', sprintf('Row %d (AUC)', rowIdx));
end
xlabel('Sheet Index');
ylabel('Estimated AUC Value');
title('Estimated AUC values for First Three Sheets');




%% Nested function to convert cells to double and non-convertible cells to NaN
function num = convertToNumeric(cellContent)
    if isnumeric(cellContent)
        num = cellContent;
    elseif isnan(str2double(cellContent))
        num = NaN;
    else
        num = str2double(cellContent);
    end
end
% Define a function for inline if-else
function result = iff(condition, true_value, false_value)
    if condition
        result = true_value;
    else
        result = false_value;
    end
end

function x_val = findXForYHalf(params)
    func = @(x) params(1) / (1 + params(2) * x) - 0.5;
    if func(0) * func(9) > 0
        % The function doesn't cross y = 0.5 within [0, 9]
        x_val = NaN;
    else
        x_val = fzero(func, [0, 9]);
    end
end