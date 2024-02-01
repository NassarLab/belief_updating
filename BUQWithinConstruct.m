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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Response CDF %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot CDFs for Each Question with Key Statistics
numQuestions = size(allData, 2); % Assuming allData has participants as rows and questions as columns

% Adjust these values based on your preference and number of questions
histSubplotsPerRow = 5;
histSubplotsPerCol = 5;
totalHistSubplots = histSubplotsPerRow * histSubplotsPerCol;

numHistFigures = ceil(numQuestions / totalHistSubplots);

for fig = 1:numHistFigures
    figure('Position', [100, 400, 1700, 1000]);
    for subPlotNum = 1:totalHistSubplots
        questionNum = (fig - 1) * totalHistSubplots + subPlotNum;
        if questionNum > numQuestions
            break;
        end
        subplot(histSubplotsPerRow, histSubplotsPerCol, subPlotNum);
        histogram(allData(:, questionNum), 'Normalization', 'probability', 'DisplayName', 'Histogram');
        hold on;
        [f, x] = ecdf(allData(:, questionNum));
        plot(x, f, 'r-', 'LineWidth', 2, 'DisplayName', 'CDF');
        hold off;
        legend;
        title(['Histogram/CDF for Q' num2str(questionNum)]);
        xlabel('Response');
        ylabel('Frequency / Cumulative Probability');

        % Key statistics
        medianVal = median(allData(:, questionNum));
        percentile25 = prctile(allData(:, questionNum), 25);
        percentile75 = prctile(allData(:, questionNum), 75);
        
        % Marking statistics on the plot
        hold on;
        xline(medianVal, '--r', 'Median');
        xline(percentile25, ':g', '25th Percentile');
        xline(percentile75, ':b', '75th Percentile');
        hold off;
        
        title(['CDF for Q' num2str(questionNum)]);
        xlabel('Response');
        ylabel('Cumulative Probability');
        grid on;
    end
end

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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%% Out-Of Construct Correlation %%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Number of Components and Questions for Analysis
% numComponentsToAnalyze = 6;
% numQuestionsPerComponent = 10;
% 
% % Initialize cell array to store unique top questions for each component
% topQuestionsPerComponent = cell(numComponentsToAnalyze, 1);
% usedQuestions = []; % Track already used questions
% 
% for i = 1:numComponentsToAnalyze
%     % Get top questions based on loadings, excluding already used ones
%     loadings = abs(coeff(:, i));
%     loadings(usedQuestions) = 0; % Exclude used questions
%     [~, sortedIndices] = sort(loadings, 'descend');
%     
%     topQuestions = sortedIndices(1:numQuestionsPerComponent);
%     topQuestionsPerComponent{i} = topQuestions;
%     usedQuestions = [usedQuestions; topQuestions]; % Update used questions
% end
% 
% % Initialize matrix to store correlations
% correlationMatrix = zeros(numComponentsToAnalyze * numQuestionsPerComponent);
% 
% % Calculate correlations
% for i = 1:numComponentsToAnalyze
%     for j = 1:numComponentsToAnalyze
%         questionsI = topQuestionsPerComponent{i};
%         questionsJ = topQuestionsPerComponent{j};
%         
%         % Extract data for the selected questions
%         dataI = allData(:, questionsI);
%         dataJ = allData(:, questionsJ);
%         
%         % Calculate correlation matrix
%         idxStart = (i-1)*numQuestionsPerComponent + 1;
%         idxEnd = i*numQuestionsPerComponent;
%         idxStartJ = (j-1)*numQuestionsPerComponent + 1;
%         idxEndJ = j*numQuestionsPerComponent;
%         
%         correlationMatrix(idxStart:idxEnd, idxStartJ:idxEndJ) = corr(dataI, dataJ);
%     end
% end
% 
% % Display and visualize the correlation matrix
% disp('Correlation Matrix:');
% disp(correlationMatrix);
% 
% % Visualization
% figure;
% imagesc(correlationMatrix);
% colorbar;
% colormap(parula);
% title('Correlation Matrix of Top Questions from Top Components');
% xlabel('Questions');
% ylabel('Questions');
% axis square;
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%% Within Construct Correlation %%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Within Construct Correlation for First 6 Components with Specified Number of Questions
% numComponentsToAnalyze = 6;
% numQuestionsPerComponent = 10; % Set the number of questions per component
% 
% % Initialize cell array to store correlation matrices for the first 6 components
% firstSixComponentCorrelationMatrices = cell(numComponentsToAnalyze, 1);
% 
% for i = 1:numComponentsToAnalyze
%     % Get the absolute loadings for the current component
%     componentLoadings = abs(coeff(:, i));
% 
%     % Sort the loadings and get the indices of the top N loadings
%     [~, sortedIndices] = sort(componentLoadings, 'descend');
%     topIndices = sortedIndices(1:min(numQuestionsPerComponent, length(sortedIndices)));
%     
%     % Check if there are enough variables for a meaningful correlation matrix
%     if length(topIndices) > 1
%         % Get correlation matrix for variables with the top N loadings on the same component
%         corrMatrix = corr(allData(:, topIndices));
%         
%         % Store and print the correlation matrix for this component
%         firstSixComponentCorrelationMatrices{i} = corrMatrix;
%         fprintf('Correlation Matrix for Component %d:\n', i);
%         disp(corrMatrix);
%     else
%         firstSixComponentCorrelationMatrices{i} = NaN; % Not enough variables for a meaningful correlation matrix
%         fprintf('Insufficient Data for Component %d\n', i);
%     end
% end


% % Constants for the subplot grid
% numRows = 2;
% numCols = 3;
% 
% % Create a figure for correlation matrices of the first 6 components
% figure('Position', [100, 100, 1700, 900]);
% sgtitle('Correlation Matrices for First 6 Components');
% 
% for i = 1:numComponentsToAnalyze
%     subplot(numRows, numCols, i);
%     
%     if ~isnan(firstSixComponentCorrelationMatrices{i})
%         imagesc(firstSixComponentCorrelationMatrices{i});
%         colorbar;
%         title(['Component ', num2str(i)]);
%         xlabel('Question Number');
%         ylabel('Question Number');
%     else
%         % Display a text placeholder for components without enough variables
%         text(0.5, 0.5, 'Insufficient Data', ...
%              'HorizontalAlignment', 'center', ...
%              'VerticalAlignment', 'middle');
%         set(gca, 'XTick', [], 'YTick', []);
%         title(['Component ', num2str(i)]);
%     end
% end





% Press space key to close figure
set(gcf, 'KeyPressFcn', @(src, evt) closeFigureOnSpace(evt));

function closeFigureOnSpace(event)
    if strcmp(event.Key, 'space')
        close(gcf);
    end
end

