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

% Assuming `coeff` and `allData` are already defined from PCA
numPCs = 8; % Number of principal components to consider
numTopQuestions = 10; % Top questions per PC

% Initialize a matrix to hold the indices of the top questions for each PC
topQuestionsPC = zeros(numTopQuestions, numPCs);
topQuestionsOrder = []; % This will hold the ordered question indices as per PCs

for pc = 1:numPCs
    [~, sortedIndices] = sort(abs(coeff(:, pc)), 'descend');
    topQuestionsPC(:, pc) = sortedIndices(1:numTopQuestions);
    topQuestionsOrder = [topQuestionsOrder; sortedIndices(1:numTopQuestions)]; % Append in order
end

% Extract the data for these top questions in the PC-wise order
selectedData = zscore(allData(:, topQuestionsOrder));

% Calculate the correlation matrix for the selected questions in the order
corrMatrixSelected = corr(selectedData);

% Visualize the correlation matrix
figure;
imagesc(corrMatrixSelected);
colorbar;
title('Correlation Matrix of Top 6 Questions for Each PC, Sorted by PC');
colormap jet;
axis square;

% Set the tick labels to reflect the PC and question order
xticks(1:length(topQuestionsOrder));
yticks(1:length(topQuestionsOrder));
xticklabels(arrayfun(@(x) ['Q', num2str(x)], topQuestionsOrder, 'UniformOutput', false));
yticklabels(arrayfun(@(x) ['Q', num2str(x)], topQuestionsOrder, 'UniformOutput', false));
xtickangle(45);

% Improve readability by adding grid lines to separate PCs
hold on;
for pc = 1:numPCs-1
    line([pc*numTopQuestions+0.5, pc*numTopQuestions+0.5], ylim, 'Color', 'w', 'LineWidth', 2);
    line(xlim, [pc*numTopQuestions+0.5, pc*numTopQuestions+0.5], 'Color', 'w', 'LineWidth', 2);
end
hold off;

% Adjust the interpreter for tick labels if needed
set(gca, 'TickLabelInterpreter', 'none');


% Initialize variables to hold average correlations
avgCorrWithinPCs = zeros(1, numPCs);
avgCorrBetweenPCs = zeros(numPCs, numPCs); % This will be a symmetric matrix

% Calculate average correlation within PCs
for pc = 1:numPCs
    startIndex = (pc-1) * numTopQuestions + 1;
    endIndex = pc * numTopQuestions;
    pcBlock = corrMatrixSelected(startIndex:endIndex, startIndex:endIndex);
    
    % Exclude the diagonal (self-correlation) when calculating the average
    avgCorrWithinPCs(pc) = abs(mean(pcBlock(triu(true(size(pcBlock)), 1)), 'all'));
end

% Calculate average correlation between PCs
for pc1 = 1:numPCs
    for pc2 = pc1+1:numPCs % Only need to calculate for one half due to symmetry
        startIdx1 = (pc1-1) * numTopQuestions + 1;
        endIdx1 = pc1 * numTopQuestions;
        startIdx2 = (pc2-1) * numTopQuestions + 1;
        endIdx2 = pc2 * numTopQuestions;
        
        % Extract the correlation block between pc1 and pc2
        betweenBlock = corrMatrixSelected(startIdx1:endIdx1, startIdx2:endIdx2);
        
        % Calculate and store the average correlation
        avgCorrBetweenPCs(pc1, pc2) = mean(betweenBlock, 'all');
        avgCorrBetweenPCs(pc2, pc1) = avgCorrBetweenPCs(pc1, pc2); % Mirror due to symmetry
    end
end

% Display the results
disp('Average correlation within PCs:');
disp(avgCorrWithinPCs);


%%%%%%% Graph 2 %%%%%%%%
% Initialize variables for the absolute correlations
minAbsBetweenPCs = zeros(1, numPCs);
maxAbsBetweenPCs = zeros(1, numPCs);
meanAbsBetweenPCs = zeros(1, numPCs);

% Calculate the minimum, maximum, and mean of the absolute correlations
for pc = 1:numPCs
    betweenCorrs = abs(avgCorrBetweenPCs(pc, :)); % Take absolute values
    betweenCorrs(pc) = []; % Remove self-comparison
    minAbsBetweenPCs(pc) = min(betweenCorrs);
    maxAbsBetweenPCs(pc) = max(betweenCorrs);
    meanAbsBetweenPCs(pc) = mean(betweenCorrs);
end

% Visualization with average absolute correlations
figure;
scatter(1:numPCs, avgCorrWithinPCs, 100, 'Filled', 'DisplayName', 'Avg Within-PC Correlation', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', [0 0.5 0.5]);
hold on;

% Plot lines to show the range of absolute between-PC correlations for each PC
for pc = 1:numPCs
    line([pc, pc], [minAbsBetweenPCs(pc), maxAbsBetweenPCs(pc)], 'Color', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
end

% Add error bars to show the range of absolute between-PC correlations (from min to max)
errorbar(1:numPCs, meanAbsBetweenPCs, minAbsBetweenPCs - meanAbsBetweenPCs, maxAbsBetweenPCs - meanAbsBetweenPCs, 'k', 'linestyle', 'none', 'DisplayName', 'Absolute Between-PC Range');

% Optionally, plot the mean absolute between-PC correlations as a separate line or markers
plot(1:numPCs, meanAbsBetweenPCs, 'r-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Mean Absolute Between-PC');

hold off;
title('Comparison of Correlation Within and Between PCs (Absolute Values)');
xlabel('Principal Component (PC)');
ylabel('Average Absolute Correlation');
legend('Location', 'Best');
set(gca, 'XTick', 1:numPCs, 'XTickLabel', arrayfun(@(x) ['PC', num2str(x)], 1:numPCs, 'UniformOutput', false));

% % Parameters for customization
% numPCs = 10; % Number of principal components to analyze
% numScores = 10; % Number of top scores to select based on magnitude
% 
% % Extract scores for the specified number of PCs
% scoresPCs = score(:, 1:numPCs);
% 
% % Initialize matrix to hold indices and scores of top elements
% topIndices = zeros(numScores, numPCs); % Indices of top scores
% topScoresPCs = zeros(numScores, numPCs); % Top scores for each PC
% 
% % Select the top elements
% for i = 1:numPCs
%     [~, sortedIndices] = sort(abs(scoresPCs(:, i)), 'descend');
%     topIndices(:, i) = sortedIndices(1:numScores);
%     topScoresPCs(:, i) = scoresPCs(topIndices(:, i), i);
% end
% 
% % Calculate the correlation matrix among the top scores of the specified PCs
% corrMatrixTopScores = corr(topScoresPCs);
% 
% % Display the correlation matrix in the console
% disp(['Correlation matrix among the top ', num2str(numScores), ' scores of first ', num2str(numPCs), ' PCs:']);
% disp(corrMatrixTopScores);
% 
% % Visualize the correlation matrix as a heatmap
% figure;
% imagesc(corrMatrixTopScores); 
% colorbar;
% title(['Correlation Matrix of Top ', num2str(numScores), ' Questions for First ', num2str(numPCs), ' PCs']);
% colormap; 
% set(gca, 'XTick', 1:numPCs, 'YTick', 1:numPCs, ...
%     'XTickLabels', arrayfun(@(x) ['PC', num2str(x)], 1:numPCs, 'UniformOutput', false), ...
%     'YTickLabels', arrayfun(@(x) ['PC', num2str(x)], 1:numPCs, 'UniformOutput', false), ...
%     'TickLabelInterpreter', 'none');
% axis square; 
% 
% % Adjust the number of PCs based on the length of averageInFactorCorrelation
% numPCs = length(averageInFactorCorrelation);  % This ensures the loop does not exceed the array bounds
% 
% % Graph
% figure;
% hold on; % Hold on to plot multiple data points on the same graph
% 
% % Plot settings
% markers = {'s', '^'}; % Square and triangle markers
% colors = {'b', 'r'}; % Blue for in-factor, red for out-factor
% 
% % Iterate only within the bounds of averageInFactorCorrelation
% for i = 1:numPCs
%     % Plot in-factor correlation coefficient with square marker
%     plot(i, dataToPlot(1, i), markers{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'MarkerSize', 8);
%     
%     % Plot out-factor correlation coefficient with triangle marker
%     plot(i, dataToPlot(2, i), markers{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'MarkerSize', 8);
% end
% 
% % Customize the graph
% set(gca, 'XTick', 1:numPCs, 'XTickLabel', arrayfun(@(x) ['PC', num2str(x)], 1:numPCs, 'UniformOutput', false));
% ylabel('Average Correlation Coefficient');
% title('In-Factor and Out-Factor Average Correlation Coefficients for PCs');
% legend({'In-Factor', 'Out-Factor'}, 'Location', 'BestOutside');
% axis square; 
% hold off; 
% 
% 
% % Customize the graph
% set(gca, 'XTick', 1:numPCs, 'XTickLabel', arrayfun(@(x) ['PC', num2str(x)], 1:numPCs, 'UniformOutput', false));
% ylabel('Average Correlation Coefficient');
% title('In-Factor and Out-Factor Average Correlation Coefficients for PCs');
% legend({'In-Factor', 'Out-Factor'}, 'Location', 'BestOutside');
% axis square; 
% hold off; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% With-in Construct Correlation %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Assume 'score' is the matrix containing your PCA scores
% % Extract the scores for the first 6 principal components
% PCs = score(:, 1:6);
% 
% % Calculate the correlation matrix among these six principal components
% corrMatrix = corr(PCs);
% 
% % Display the correlation matrix
% disp('Correlation matrix among PC1, PC2, PC3, PC4, PC5, and PC6:');
% disp(corrMatrix);
% 
% % Optionally, visualize the correlation matrix as a heatmap
% figure;
% imagesc(corrMatrix);
% colorbar;
% title('Correlation Matrix among PC1 to PC6');
% xticks(1:6);
% yticks(1:6);
% xticklabels({'PC1','PC2','PC3','PC4','PC5','PC6'});
% yticklabels({'PC1','PC2','PC3','PC4','PC5','PC6'});
% axis square;




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

