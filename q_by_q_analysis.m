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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% PCA %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% Choose the principal component to plot, e.g., the first principal component
pc1Loadings = coeff(:, 1); 

% Create a figure for the loadings of the chosen principal component
figure;
bar(pc1Loadings);
title('Loadings of Questions on the First Principal Component (PC1)');
xlabel('Question');
ylabel('Loading Value');

% Adjust the x-axis to show the full, cleaned question texts as labels
% set(gca, 'XTick', 1:length(cleanedQuestions), 'XTickLabel', cleanedQuestions);
xtickangle(45); % Rotate labels for better readability
set(gca, 'FontSize', 8); % Adjust font size for visibility
grid on;

% Optionally, highlight significant loadings on the graph
% Let's assume significant loadings are those with absolute values above a certain threshold, e.g., 0.2
hold on;
significantIdx = find(abs(pc1Loadings) > 0.2);
bar(significantIdx, pc1Loadings(significantIdx), 'r');
legend('Loadings', 'Significant Loadings (>0.2)');
hold off;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Response CDF %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot CDFs for Each Question with Key Statistics
% numQuestions = size(allData, 2); % Assuming allData has participants as rows and questions as columns

% % Adjust these values based on your preference and number of questions
% histSubplotsPerRow = 5;
% histSubplotsPerCol = 5;
% totalHistSubplots = histSubplotsPerRow * histSubplotsPerCol;
% 
% numHistFigures = ceil(numQuestions / totalHistSubplots);
% 
% for fig = 1:numHistFigures
%     figure('Position', [100, 400, 1700, 1000]);
%     for subPlotNum = 1:totalHistSubplots
%         questionNum = (fig - 1) * totalHistSubplots + subPlotNum;
%         if questionNum > numQuestions
%             break;
%         end
%         subplot(histSubplotsPerRow, histSubplotsPerCol, subPlotNum);
%         histogram(allData(:, questionNum), 'Normalization', 'probability', 'DisplayName', 'Histogram');
%         hold on;
%         [f, x] = ecdf(allData(:, questionNum));
%         plot(x, f, 'r-', 'LineWidth', 2, 'DisplayName', 'CDF');
%         hold off;
%         legend;
%         title(['Histogram/CDF for Q' num2str(questionNum)]);
%         xlabel('Response');
%         ylabel('Frequency / Cumulative Probability');
% 
%         % Key statistics
%         medianVal = median(allData(:, questionNum));
%         percentile25 = prctile(allData(:, questionNum), 25);
%         percentile75 = prctile(allData(:, questionNum), 75);
%         
%         % Marking statistics on the plot
%         hold on;
%         xline(medianVal, '--r', 'Median');
%         xline(percentile25, ':g', '25th Percentile');
%         xline(percentile75, ':b', '75th Percentile');
%         hold off;
%         
%         title(['CDF for Q' num2str(questionNum)]);
%         xlabel('Response');
%         ylabel('Cumulative Probability');
%         grid on;
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Skewness %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Key Statistics for Each Question
numQuestions = size(allData, 2);
medians = zeros(numQuestions, 1);
percentile25s = zeros(numQuestions, 1);
percentile75s = zeros(numQuestions, 1);
minVals = zeros(numQuestions, 1);
maxVals = zeros(numQuestions, 1);
skewnessVals = zeros(numQuestions, 1);

for i = 1:numQuestions
    medians(i) = median(allData(:, i));
    percentile25s(i) = prctile(allData(:, i), 25);
    percentile75s(i) = prctile(allData(:, i), 75);
    minVals(i) = min(allData(:, i));
    maxVals(i) = max(allData(:, i));
    skewnessVals(i) = skewness(allData(:, i));
end

questions = data.Question;

%% Print Questions with Skewness Higher than 0.5 or Lower than -0.5
highSkewnessQuestions = find(skewnessVals > 0.5);
lowSkewnessQuestions = find(skewnessVals < -0.5);

fprintf('Questions with skewness higher than 0.5:\n');
for i = 1:length(highSkewnessQuestions)
    qNum = highSkewnessQuestions(i);
    fprintf('Question %d (%.2f): %s\n', qNum, skewnessVals(qNum), questions{qNum});
end

fprintf('\nQuestions with skewness lower than -0.5:\n');
for i = 1:length(lowSkewnessQuestions)
    qNum = lowSkewnessQuestions(i);
    fprintf('Question %d (%.2f): %s\n', qNum, skewnessVals(qNum), questions{qNum});
end

%% Modify Question Texts for Skewness Plot
% Clean "<b>" from the beginning of each question text
cleanedQuestions = strrep(questions, '<b>', '');

% Generate a skewness plot
figure; % Creates a new figure
bar(skewnessVals); % Creates a bar graph of the skewness values
title('Skewness of Responses for Each Question');
xlabel('Question Number');
ylabel('Skewness Value');
xticks(1:length(cleanedQuestions)); % Adjusts the x-axis ticks to match the number of questions
xtickangle(45); % Rotates the x-axis labels for better readability
xticklabels(cleanedQuestions); % Use cleaned questions as x-axis labels
grid on; % Adds a grid for easier visualization of skewness values

% Enhance the plot to highlight high and low skewness
hold on;
highSkewnessIdx = find(skewnessVals > 0.5);
lowSkewnessIdx = find(skewnessVals < -0.5);
plot(highSkewnessIdx, skewnessVals(highSkewnessIdx), 'r*', 'MarkerSize', 10, 'LineWidth', 2);
plot(lowSkewnessIdx, skewnessVals(lowSkewnessIdx), 'b*', 'MarkerSize', 10, 'LineWidth', 2);
legend('Skewness', 'High Skewness (>0.5)', 'Low Skewness (<-0.5)');
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Variance, Mean, Median %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize matrices to hold statistics
numQuestions = size(allData, 2);
variances = zeros(numQuestions, 1);
means = zeros(numQuestions, 1);
medians = zeros(numQuestions, 1);
responseCounts = zeros(numQuestions, 1);

% Calculate statistics for each question
for i = 1:numQuestions
    currentData = allData(:, i); % Extract data for current question
    nonMissingData = currentData(~isnan(currentData)); % Exclude missing responses if any
    variances(i) = var(nonMissingData);
    means(i) = mean(nonMissingData);
    medians(i) = median(nonMissingData);
    responseCounts(i) = length(nonMissingData); % Count of non-missing responses
end

% Rank the statistics
[~, rankVariances] = sort(variances, 'descend');
[~, rankMeans] = sort(means, 'descend');
[~, rankMedians] = sort(medians, 'descend');
[~, rankResponseCounts] = sort(responseCounts, 'descend');

% Create a table to display the ranked questions for each statistic
rankedTable = table(rankVariances, rankMeans, rankMedians, rankResponseCounts, ...
    'VariableNames', {'Variance', 'Mean', 'Median', 'ResponseCount'});

% Optionally, display the table
disp(rankedTable);

% To access rankings for a specific question, e.g., question number 120
qNum = 120;
fprintf('Rankings for Question %d: Variance: %d, Mean: %d, Median: %d, Response Count: %d\n', ...
    qNum, find(rankVariances == qNum), find(rankMeans == qNum), find(rankMedians == qNum), find(rankResponseCounts == qNum));

% Assuming 'questions' variable contains question texts
cleanedQuestions = strrep(questions, '<b>', ''); % Remove "<b>" from questions
[sortedVariances, sortedIdxVariances] = sort(variances, 'descend');

%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Graph %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
% figure;
% plot(sort(variances, 'descend'));
% hold on;
% plot(means(sortedIdxVariances));
% 
% title('Variance of Responses for Each Question');
% xlabel('Question Number');
% ylabel('Variance');
% set(gca, 'XTick', 1:length(cleanedQuestions), 'XTickLabel', cleanedQuestions); % Use cleaned questions as x-axis labels
% xtickangle(45); % Rotate labels for better readability
% set(gca, 'FontSize', 8); % Adjust font size here
% grid on;
% 
% % Highlighting the top N questions with highest variance
% [sortedVariances, sortedIdxVariances] = sort(variances, 'descend');
% hold on;
% bar(sortedIdxVariances(1:5), sortedVariances(1:5), 'r');
% legend('All Questions', 'Top 5 Variances');
% hold off;
% 
% 
% % Highlighting the top N questions with highest variance
% [sortedVariances, sortedIdxVariances] = sort(variances, 'descend');
% hold on;
% bar(sortedIdxVariances(1:5), sortedVariances(1:5), 'r');
% legend('All Questions', 'Top 5 Variances');
% hold off;
% 
% % Mean
% figure;
% bar(means);
% title('Mean of Responses for Each Question');
% xlabel('Question Number');
% ylabel('Mean');
% set(gca, 'XTick', 1:length(cleanedQuestions), 'XTickLabel', cleanedQuestions);
% xtickangle(45);
% grid on;
% 
% % Highlighting the top N questions with highest mean
% [sortedMeans, sortedIdxMeans] = sort(means, 'descend');
% hold on;
% bar(sortedIdxMeans(1:5), sortedMeans(1:5), 'r');
% legend('All Questions', 'Top 5 Means');
% hold off;
% 
% % Median
% figure;
% bar(medians);
% title('Median of Responses for Each Question');
% xlabel('Question Number');
% ylabel('Median');
% set(gca, 'XTick', 1:length(cleanedQuestions), 'XTickLabel', cleanedQuestions);
% xtickangle(45);
% grid on;
% 
% % Highlighting the top N questions with highest median
% [sortedMedians, sortedIdxMedians] = sort(medians, 'descend');
% hold on;
% bar(sortedIdxMedians(1:5), sortedMedians(1:5), 'r');
% legend('All Questions', 'Top 5 Medians');
% hold off;
% 
% 
% % Response count
% figure;
% bar(responseCounts);
% title('Number of Responses for Each Question');
% xlabel('Question Number');
% ylabel('Response Count');
% set(gca, 'XTick', 1:length(cleanedQuestions), 'XTickLabel', cleanedQuestions);
% xtickangle(45);
% grid on;
% 
% % Highlighting the questions with highest response counts
% [sortedCounts, sortedIdxCounts] = sort(responseCounts, 'descend');
% hold on;
% bar(sortedIdxCounts(1:5), sortedCounts(1:5), 'r');
% legend('All Questions', 'Top 5 Response Counts');
% hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% bottom 20% %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assuming variances, means, and medians are already calculated as shown above
numQuestions = size(allData, 2); % Total number of questions
bottomPercentCount = ceil(numQuestions * 0.5); % Number of questions in the bottom %

% Sort variances to find the bottom 10%
[sortedVariances, sortIndex] = sort(variances, 'ascend');
bottomPercentIndex = sortIndex(1:bottomPercentCount);

% Extract statistics for the bottom 10%
bottomVariances = variances(bottomPercentIndex);
bottomMeans = means(bottomPercentIndex);
bottomMedians = medians(bottomPercentIndex);

% Plotting
figure;
hold on;
plot(1:bottomPercentCount, bottomVariances, 'b-o', 'LineWidth', 2, 'DisplayName', 'Variance');
plot(1:bottomPercentCount, bottomMeans, 'r-*', 'LineWidth', 2, 'DisplayName', 'Mean');
plot(1:bottomPercentCount, bottomMedians, 'g-s', 'LineWidth', 2, 'DisplayName', 'Median');
hold off;

% Enhance the plot
title('Bottom 20% Questions by Variance: Variance, Mean, and Median Responses');
xlabel('Question Index');
ylabel('Value');
legend('show');
grid on;

% Adjust the x-ticks to reflect the actual question numbers
xticks(1:bottomPercentCount);
xticklabels(arrayfun(@(x) sprintf('Q%d', x), bottomPercentIndex, 'UniformOutput', false));
xtickangle(45);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% by construct %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define question numbers for each construct
constructs = struct(...
    'Nuance', [17, 27, 47, 79, 94, 101, 111, 124, 132, 141], ...
    'FalseBeliefs', [7, 12, 119, 23, 34, 40, 49, 73, 81, 86, 91, 98, 101, 103, 109, 120, 124, 128, 136, 139, 142, 145, 147], ...
    'Evidence', [1, 16, 22, 26, 33, 37, 45, 52, 54, 58, 64, 72, 76, 79, 82, 85, 89, 94, 97, 100, 104, 106, 111, 117, 119, 125, 131, 133, 137, 141, 144, 150], ...
    'ConflictingBeliefs', [3, 21, 32, 36, 43, 57, 62, 68, 87, 108, 123], ...
    'Politics', [4, 6, 9, 13, 27, 29, 38, 41, 45, 50, 55, 64, 65, 69, 76, 81, 91, 108, 111, 113, 115, 119, 126, 146], ...
    'MoralsAndTruth', [2, 10, 49, 60, 68, 75, 97], ...
    'ReasoningStyles', [8, 14, 24, 30, 37, 43, 54, 59, 62, 73]);
% 'AttentionCheck', [22, 71, 115], ...
constructNames = fieldnames(constructs);
constructStats = struct();

% Loop through each construct
% Loop through each construct
% Loop through each construct
for c = 1:length(constructNames)
    construct = constructNames{c};
    questionIndices = constructs.(construct);

    % Ensure questionIndices are sorted in ascending order for processing
    questionIndices = sort(questionIndices(questionIndices <= size(allData, 2)));

    % Skip if no valid questions
    if isempty(questionIndices)
        continue;
    end

    % Initialize arrays to store statistics for each question in the construct
    questionVariances = zeros(1, length(questionIndices));
    questionMeans = zeros(1, length(questionIndices));
    questionMedians = zeros(1, length(questionIndices));

    % Calculate statistics for each question
    for q = 1:length(questionIndices)
        qIndex = questionIndices(q);
        questionData = allData(:, qIndex);

        questionVariances(q) = var(questionData, 'omitnan');
        questionMeans(q) = mean(questionData, 'omitnan');
        questionMedians(q) = median(questionData, 'omitnan');
    end

    % Sort the questions by variance in descending order
    [sortedVariances, sortIndices] = sort(questionVariances, 'descend');
    sortedQuestionIndices = questionIndices(sortIndices);
    sortedQuestionMeans = questionMeans(sortIndices);
    sortedQuestionMedians = questionMedians(sortIndices);

    % Plotting the statistics for the current construct in sorted order
    figure; % Create a new figure for each construct
    hold on;

    % Plot each statistic using the sorted order
    plot(1:length(sortedQuestionIndices), sortedVariances, 'b-o', 'LineWidth', 2, 'DisplayName', 'Variance');
    plot(1:length(sortedQuestionIndices), sortedQuestionMeans, 'r-*', 'LineWidth', 2, 'DisplayName', 'Mean');
    plot(1:length(sortedQuestionIndices), sortedQuestionMedians, 'g-s', 'LineWidth', 2, 'DisplayName', 'Median');

    hold off;
    title(['Questions in Construct: ', construct]);
    xlabel('Sorted Question Index');
    ylabel('Statistic Value');
    legend('show');
    grid on;

    % Adjust the x-ticks to reflect the sorted question indices and use custom labels
    xticks(1:length(sortedQuestionIndices));
    xticklabels(arrayfun(@(x) strcat('Q', num2str(x)), sortedQuestionIndices, 'UniformOutput', false));
    xtickangle(45);
end








% Press space key to close figure
set(gcf, 'KeyPressFcn', @(src, evt) closeFigureOnSpace(evt));

function closeFigureOnSpace(event)
    if strcmp(event.Key, 'space')
        close(gcf);
    end
end


