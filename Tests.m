% Clear workspace and window
clear; clc;

% Choose if you would like to see the results displayed in a single figure 
% window ('subplots') or in multiple figure windows ('manyplots')
dispres = 'manyplots';

LW = 2;  % Plot line width
MS = 10; % Size of markers on plots

numTests = 5;
SNR = 10.^(0:8)';

model = 'CT';  % Choose between 'blur' and 'CT'
n_values = [64, 128];  % List of n values to iterate over

resultsDir = fullfile(fileparts(mfilename('fullpath')), 'Results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

for n = n_values
  % Define the test problem
  [A, b, x, ProbInfo, m0] = defineTestProblem(model, n);
  % Plot singular values decay of A
  figure(5); clf;
  [~,S,~] = svds(A, m0);
  loglog(diag(S), 'linewidth', LW);
  title('Singular values decay of A', 'interpreter', 'latex', 'fontsize', 18);

  NoiseLevel = (norm(b) / sqrt(size(b,1))) ./ SNR;
  % Preallocate table for test data
  testData = table('Size', [0, 3], 'VariableTypes', {'double', 'double', 'string'},...
                   'VariableNames', {'SNR', 'Error', 'Method'});

  for i = 1:size(SNR,1)
    for j = 1:numTests
      fprintf('Test %d, SNR = %f\n', j, SNR(i));
      rng(j);  % Set seed for reproducibility

      [bn, NoiseInfo] = PRnoise(b, 'gauss', NoiseLevel(i));
      [X_gcv, X_opt, error_gcv, error_opt] = gcv(A, x, bn, m0);
      testData = [testData; table(SNR(i), error_gcv, "gcv",...
                                  'VariableNames', {'SNR', 'Error', 'Method'})];
      testData = [testData; table(SNR(i), error_opt, "opt",...
                                  'VariableNames', {'SNR', 'Error', 'Method'})];
    end
  end
  % Plot results
  plotResults(testData, SNR, n);

  % Display the reconstructions
  displayReconstructions(x, b, X_gcv, X_opt, ProbInfo);

  % Save figures if needed
  saveFigures(dispres, n, resultsDir, model);
end



% Function to define the test problem
function [A, b, x, ProbInfo, m0] = defineTestProblem(model, n)
  if strcmp(model, 'blur')
    m0 = n*n;
    [A, b, x, ProbInfo] = PRblurrotation(n);
  elseif strcmp(model, 'CT')
    options = IRset();
    options.sm = true;
    if n == 64
      m0 = 4096;
    elseif n == 128
      m0 = 16384;
    end
    [A, b, x, ProbInfo] = PRtomo(n, options);
  end
end

% Function to plot results
function plotResults(testData, SNR, n)
  figure(10); clf;
  testData.SNR = categorical(testData.SNR, SNR,...
   {'1', '1e1', '1e2', '1e3', '1e4', '1e5', '1e6', '1e7', '1e8'});
  boxchart(testData.SNR, testData.Error, 'GroupByColor', testData.Method);
  title(['Results for n = ', num2str(n)]);
  ylabel('Relative error');
  xlabel('SNR');
  ax = gca;
  ax.YAxis.Scale ="log";
  legend()
end

% Function to display reconstructions
function displayReconstructions(x, b, X_gcv, X_opt, ProbInfo)
  figure(1); clf;
  PRshowx(x, ProbInfo);
  set(gca, 'fontsize', 24);
  title('True solution', 'interpreter', 'latex', 'fontsize', 18);

  figure(2); clf;
  PRshowb(b, ProbInfo);
  set(gca, 'fontsize', 24);
  title('Noisy data', 'interpreter', 'latex', 'fontsize', 18);

  figure(3); clf;
  PRshowx(X_gcv, ProbInfo);
  title('GCV sol.', 'interpreter', 'latex', 'fontsize', 18);

  figure(4); clf;
  PRshowx(X_opt, ProbInfo);
  title('Optimal sol.', 'interpreter', 'latex', 'fontsize', 18);
end

% Function to save figures
function saveFigures(dispres, n, resultsDir, model)
  if strcmp(dispres, 'subplots')
    saveFigure(fullfile(resultsDir, ['AllPlots_n', num2str(n), model, '.png']), 1);
  elseif strcmp(dispres, 'manyplots')
    saveFigure(fullfile(resultsDir, ['Orig_n', num2str(n), model, '.eps']), 1);
    saveFigure(fullfile(resultsDir, ['Signal_n', num2str(n), model, '.eps']), 2);
    saveFigure(fullfile(resultsDir, ['GCV_n', num2str(n), model, '.eps']), 3);
    saveFigure(fullfile(resultsDir, ['Opt_n', num2str(n), model, '.eps']), 4);
    saveFigure(fullfile(resultsDir, ['SVdecay_n', num2str(n), model, '.pdf']), 5);
    saveFigure(fullfile(resultsDir, ['Errors_n', num2str(n), model, '.pdf']), 10);
  end
end

% Helper function to save a figure
function saveFigure(filename, figNumber)
  exportgraphics(figure(figNumber), filename);
end