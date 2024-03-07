%% ======== Belief Updating Questionnaire PCA Analysis script

%% Read in Data
% Set the directory path
directoryPath = 'buq_response_filtered';
% Bad data, took these out manually 
% 1379, 1416, 1422, 1433, 1437, 1443. 1449, 1458, 1476, 1479,
% 1500, 1503, 1506, 1513(no 1), 1515, 1518, 1524, 1528, 1529, 

% Open file to write
resultsFile = fopen('PCA_Results.txt', 'w');
if resultsFile == -1
    error('Cannot open PCA_Results.txt for writing.');
end

% Get a list of all .csv files in the directory
csvFiles = dir(fullfile(directoryPath, '*.csv'));

% Initialize a cell array to store all tables
allTables = cell(length(csvFiles), 1);

% Initialize an empty struct to store all tables
% Loop through each file and read it as a table
for i = 1:length(csvFiles)
    fullFileName = fullfile(directoryPath, csvFiles(i).name);
    
    % Read the .csv file as a table
    tbl = readtable(fullFileName);
    allTables{i} = tbl;
    
    % Get a valid struct field name based on the file name
    fieldName = matlab.lang.makeValidName(csvFiles(i).name, 'ReplacementStyle', 'delete');
    
    % Store the table in the struct using the field name
    allTablesStruct.(fieldName) = tbl;
  
end


%% PCA
a = fieldnames(allTablesStruct);

allData=[];
% Initialize a logical array to store attention check results
attentionCheckResults = false(length(a), 1);  

% Parameters for non-engagement detection
% variance threshold
varianceThreshold = 0.5;  
sequentialSameAnswerThreshold = 8;
% Threshold for straight lining detection, e.g., 90% same answers
straightLineThreshold = 0.8;

for i = 1:length(a)
    fprintf('processing: %s\n', a{i});
    
    try  
        eval(sprintf('data=allTablesStruct.%s;', a{i}));

        if (length(data.answer_num) ~= 150)
            fprintf('Format Error, check %s', a{i})
            continue;
        end
    
        % Check attention check questions
        attentionCheckResults(i) = data.answer_num(22) == 3 && ...
                               data.answer_num(71) == 1 && ...
                               data.answer_num(115) == 1;
    
        % Attention check filter
        if attentionCheckResults(i)
            tempData = data.answer_num;
            tempData([22, 71, 115]) = []; % Remove attention check questions
            
            
            % Sequential same answer check
            sequentialSameAnswerCount = max(diff([0; find(diff(tempData)); numel(tempData)]));
            if sequentialSameAnswerCount > sequentialSameAnswerThreshold
                fprintf('Sequential same answers detected for %s\n', a{i});
                fprintf('%s skipped\n \n', a{i});
                continue;  % Skip to the next iteration
            end
            
            % Straight-lining check
            modeAnswer = mode(tempData);
            if mean(tempData == modeAnswer) > straightLineThreshold
                fprintf('Straight-lining detected for %s\n', a{i});
                continue;  % Skip to the next iteration
            end

             % Variance check
            if var(tempData) < varianceThreshold
                fprintf('Low variance detected for %s\n', a{i});
                %continue;
            end
       
            % Write data if pass all checks
            allData(end+1, :)= tempData;

        else
            fprintf('Attention check failed for %s: Q22=%d, Q71=%d, Q115=%d\n', ...
                a{i}, data.answer_num(22), data.answer_num(71), data.answer_num(115));
        end

    catch ME
        fprintf('Error processing participant %s: %s\n', a{i}, ME.message);
    end 
end

if isempty(allData)
    error('allData is empty. No participant data passed the filters.');
end


%% Analysis
% Perform PCA on the correlation matrix
[coeff, score, latent, ~, VarExplained] = pca(zscore(allData), 'Centered', false);

% Calculate and display the correlation matrix
correlationMatrix = corr(allData);

% Display the top 'x' number of component
numComponentsToShow = 6;  
[numQuestions, numComponents] = size(coeff);

% Update the questions cell array to remove the attention check questions
questions = data.Question;

questions([22, 71, 115]) = [];

%% Writing out Result
% Displaying questions
for i = 1:numComponentsToShow
    [sortedLoadings, idx] = sort(abs(coeff(:, i)), 'descend');
    fprintf('Principal Component %d:\n', i);
    for j = 1:5  % Show top 5 influential questions for each component
        fprintf('Question %d (%.3f): %s\n', idx(j), coeff(idx(j), i), questions{idx(j)});
    end
    fprintf('\n');
end

% % Uncomment to write results to file
% for i = 1:numComponentsToShow
%     [sortedLoadings, idx] = sort(abs(coeff(:, i)), 'descend');
%     fprintf(resultsFile, 'Principal Component %d:\n', i);
%     for j = 1:5  % Show top 5 influential questions for each component
%         fprintf(resultsFile, 'Question %d (%.3f): %s\n', idx(j), coeff(idx(j), i), questions{idx(j)});
%     end
%     fprintf(resultsFile, '\n');
% end

% calculate correlation
% [r,p]=corr(allData(:,144), allData(:,104))
% [r,p]=corr(allData(:,[144, 104, 73, 10, 44]), allData(:,104))

% Close files
fclose(resultsFile);


% Define figure properties
numFigures = 5;  % Number of figures you want to display
figWidth = 400;  % Width 
figHeight = 350; % Height 
padding = 10;    % Gap between figures

% Calculate total width of all figures including padding
totalWidth = numFigures * figWidth + (numFigures - 1) * padding;

% Get the screen size
screenSize = get(0, 'ScreenSize');
screenWidth = screenSize(3);
screenHeight = screenSize(4);

% Calculate the starting x position
startX = (screenWidth - totalWidth) / 2;
startY = (screenHeight - figHeight) / 2; % You can adjust this as needed

% Define figure positions
positions = zeros(numFigures, 4);
for i = 1:numFigures
    positions(i,:) = [startX + (i-1) * (figWidth + padding), startY, figWidth, figHeight];
end


%% Figures
% Create, display, and save figures
figureFileNames = {'ParticipantResponses.png', 'CorrelationMatrix.png', 'AnotherFileName.png', 'YetAnotherFileName.png'};
 
for i = 1:numFigures
    fig = figure('Position', positions(i,:)); % Create figure with specified position
    switch i
        case 1
            imagesc(allData);
            ylabel('Participant');
            xlabel('Question');
            title('Participant Responses Matrix');
        case 2
            imagesc(correlationMatrix, [-1, 1]);  % Ensure the color scale is between -1 and 1
            title('Question Correlation Matrix');
            xlabel('Question');
            ylabel('Question');
            colorbar;
        case 3
            plot(mean(allData));
            xlabel('question');
            ylabel('average response');
            title('Question Average Response');
        case 4
            % Box Plots for each of the first six principal components on the same plane
            scoresToPlot = score(:, 1:numComponentsToShow);
            boxplot(scoresToPlot, 'Labels', arrayfun(@(x) sprintf('PC %d', x), ...
                1:numComponentsToShow, 'UniformOutput', false));
            title('Box Plots for the First Six Principal Components');
        case 5
            % Display PCA results
            subplot(2, 1, 1);
            plot(VarExplained, 'o-');
            xlabel('Principal Component');
            ylabel('Variance Explained');
            title('Scree Plot');

            subplot(2, 1, 2);
            imagesc(coeff);
            colorbar;
            xlabel('Principal Component');
            ylabel('Question');
            title('PCA Loadings');
    end
    % Save figure to file
%     saveas(fig, figureFileNames{i});
end


%% Press space key to close all figures
set(gcf, 'KeyPressFcn', @(src, evt) closeFiguresOnSpace(evt));

figures = findobj('Type', 'figure');
for i = 1:length(figures)
    set(figures(i), 'KeyPressFcn', @(src, evt) closeFiguresOnSpace(evt));
end

function closeFiguresOnSpace(event)
    if strcmp(event.Key, 'space')
        close all;
    end
end






% % Uncomment for boring printing
% figure(1)
% imagesc(allData)
% ylabel('Participant')
% xlabel('Question')
% title('Participant Responses Matrix');
% 
% figure(2)
% plot(mean(allData))
% xlabel('question')
% ylabel('average response')
% title('Question Average Response')
% 
% figure(3);
% imagesc(correlationMatrix, [-1, 1]);  % Ensure the color scale is between -1 and 1
% title('Question Correlation Matrix');
% xlabel('Question');
% ylabel('Question');
% colorbar;
% 
% % Box Plots for each of the first six principal components on the same plane
% scoresToPlot = score(:, 1:numComponentsToShow);
% 
% figure(4); % Create a single figure for the box plot
% boxplot(scoresToPlot, 'Labels', arrayfun(@(x) sprintf('PC %d', x), 1:numComponentsToShow, 'UniformOutput', false));
% title('Box Plots for the First Six Principal Components');
% 
% % Display PCA results
% figure(5);
% subplot(2, 1, 1);
% plot(VarExplained, 'o-');
% xlabel('Principal Component');
% ylabel('Variance Explained');
% title('Scree Plot');
% 
% subplot(2, 1, 2);
% imagesc(coeff);
% colorbar;
% xlabel('Principal Component');
% ylabel('Question');
% title('PCA Loadings');







% Now lets look at which questions correlate with one another:
% figure(2)
% [r]=corr(allData)
% imagesc(r, [-1, 1])
% ylabel('Question')
% xlabel('Question')
% colorbar
