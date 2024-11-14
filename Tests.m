% Clear workspace and window
clear; clc;

% Choose if you would like to see the results displayed in a single figure 
% window ('subplots') or in multiple figure windows ('manyplots')
dispres = 'manyplots';

LW = 2;  % Plot line width
MS = 10; % Size of markers on plots

numTests = 5;
SNR = 10.^(0:4)';

model = 'blur';  % Choose between 'blur' and 'CT'
n = 32;

% Preallocate table for test data
testData = table('Size', [0, 3], 'VariableTypes', {'double', 'double', 'string'},...
                 'VariableNames', {'SNR', 'Error', 'Method'});

% Define the test problem
[A, b, x, ProbInfo, m0] = defineTestProblem(model, n);

NoiseLevel = (norm(b) / sqrt(size(b,1))) ./ SNR;

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
plotResults(testData, SNR);

% Display the reconstructions
displayReconstructions(x, b, X_gcv, X_opt, ProbInfo);

% Save figures if needed
saveFigures(dispres);

% Function to define the test problem
function [A, b, x, ProbInfo, m0] = defineTestProblem(model, n)
  if strcmp(model, 'blur')
    m0 = 100;
    [A, b, x, ProbInfo] = PRblurrotation(n);
  elseif strcmp(model, 'CT')
    options = IRset();
    options.sm = true;
    m0 = 200;
    [A, b, x, ProbInfo] = PRtomo(n, options);
  end
end

% Function to plot results
function plotResults(testData, SNR)
  figure(10); clf;
  testData.SNR = categorical(testData.SNR, SNR, {'1', '1e1', '1e2', '1e3', '1e4'});
  boxchart(testData.SNR, testData.Error, 'GroupByColor', testData.Method);
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
function saveFigures(dispres)
  currentFolder = fileparts(mfilename('fullpath'));
  cd(currentFolder);
  oldcd = cd;

  try
    cd('Results');
  catch
    mkdir('Results');
    cd('Results');
  end

  if strcmp(dispres, 'subplots')
    saveFigure('AllPlots.eps', 1);
  elseif strcmp(dispres, 'manyplots')
    saveFigure('Orig.eps', 1);
    saveFigure('Signal.eps', 2);
    saveFigure('GCV.eps', 3);
    saveFigure('Opt.eps', 4);
  end
  cd(oldcd);
end

% Helper function to save a figure
function saveFigure(filename, figNumber)
  exportgraphics(figure(figNumber), filename, 'ContentType',...
                 'vector', 'Resolution', 300);
end