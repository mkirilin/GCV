% Clear workspace and window
clear; clc;

% Choose if you would like to see the results displayed in a single figure 
% window ('subplots') or in multiple figure windows ('manyplots')
dispres = 'manyplots';

LW = 2;  % Plot line width
MS = 10; % Size of markers on plots

numTests = 10;
SNR = 10.^(0:8)';

model = 'blurGauss';  % Choose between 'blur', 'blurGauss' and 'CT'
n_values = [64];  % List of n values to iterate over

resultsDir = fullfile(fileparts(mfilename('fullpath')), 'Results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

for n = n_values
  % Define the test problem
  [A, b, x, ProbInfo, m0] = defineTestProblem(model, n);
  m = size(A,1);
  % Find all non-zero singular vectors of sparce matrix A
  m_sv = min(size(A));
  fprintf('dims of A: %d x %d\n', size(A,1), size(A,2));
  assert(m0 <= m_sv, 'm0 must be less than or equal to the number of rows of A');
  fprintf('m = %d, m0 = %d\n', m, m0);
  [U,S,V] = svds(A, m_sv);

  % Plot singular values decay of A
  figure(5); clf;
  loglog(diag(S), 'linewidth', LW);
  title('Singular values decay of A', 'interpreter', 'latex', 'fontsize', 18);

  NoiseLevel = (norm(b) / sqrt(size(b,1))) ./ SNR;
  % Preallocate table for test data
  testData = table('Size', [0, 3], 'VariableTypes', {'double', 'double', 'string'},...
                   'VariableNames', {'SNR', 'Error', 'Method'});
  kData = table('Size', [0, 3], 'VariableTypes', {'double', 'double', 'string'},...
                'VariableNames', {'SNR', 'k', 'Method'});

  for i = 1:size(SNR,1)
    for j = 1:numTests
      fprintf('Test %d, SNR = %f\n', j, SNR(i));
      rng(j);  % Set seed for reproducibility

      [bn, NoiseInfo] = PRnoise(b, 'gauss', NoiseLevel(i));
      [X_gcv, X_opt, error_gcv, error_opt, k_gcv, k_opt] =...
        gcv(U(:,1:m0), S(1:m0, 1:m0), V(:,1:m0), x, bn, m0, m);
      testData = [testData; table(SNR(i), error_gcv, "gcv",...
                                  'VariableNames', {'SNR', 'Error', 'Method'})];
      testData = [testData; table(SNR(i), error_opt, "opt",...
                                  'VariableNames', {'SNR', 'Error', 'Method'})];
      kData = [kData; table(SNR(i), k_gcv, "gcv",...
                            'VariableNames', {'SNR', 'k', 'Method'})];
      kData = [kData; table(SNR(i), k_opt, "opt",...
                            'VariableNames', {'SNR', 'k', 'Method'})];
    end
  end
  % Plot results
  plotErrorResults(testData, SNR, n, m0, m);
  plotKResults(kData, SNR, n, m0, m);

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
    m0 = 500;
    [A, b, x, ProbInfo] = PRtomo(n, options);
  elseif strcmp(model, 'blurGauss')
    options = IRset();
    options.BlurLevel = 'mild';
    m0 = 4096;
    [~, b, x, ProbInfo] = PRblurgauss(n, options);
    A = psf2A(ProbInfo.psf, 'symmetric');
  else
    error('Invalid model');
  end
end

% Function to plot results
function plotErrorResults(testData, SNR, n, m0, m)
  figure(10); clf;
  testData.SNR = categorical(testData.SNR, SNR,...
   {'1', '1e1', '1e2', '1e3', '1e4', '1e5', '1e6', '1e7', '1e8'});
  boxchart(testData.SNR, testData.Error, 'GroupByColor', testData.Method);
  title(['Results for n = ', num2str(n), ', ', num2str(m0), '/', num2str(m), ' singular values'],...
        'interpreter', 'latex', 'fontsize', 18); 
  ylabel('Relative error');
  xlabel('SNR');
  set(gca,'XTickLabels', {'1', '10^{1}', '10^{2}', '10^{3}', '10^{4}', '10^{5}', '10^{6}', '10^{7}', '10^{8}'});
  ax = gca;
  ax.YAxis.Scale ="log";
  ylim(ax, ylim(ax) + [-1,1]*range(ylim(ax)).* 0.05)
  legend()
end

function plotKResults(kData, SNR, n, m0, m)
  figure(11); clf;
  kData.SNR = categorical(kData.SNR, SNR,...
   {'1', '1e1', '1e2', '1e3', '1e4', '1e5', '1e6', '1e7', '1e8'});
  boxchart(kData.SNR, kData.k, 'GroupByColor', kData.Method);
  title(['Results for n = ', num2str(n), ', ', num2str(m0), '/', num2str(m), ' singular values'],...
        'interpreter', 'latex', 'fontsize', 18); 
  ylabel('k');
  xlabel('SNR');
  set(gca,'XTickLabels', {'1', '10^{1}', '10^{2}', '10^{3}', '10^{4}', '10^{5}', '10^{6}', '10^{7}', '10^{8}'});
  ax = gca;
  ax.YAxis.Scale ="log";
  ylim(ax, ylim(ax) + [-1,1]*range(ylim(ax)).* 0.05)
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
    saveFigure(fullfile(resultsDir, ['ks_n', num2str(n), model, '.pdf']), 11);
  end
end

% Helper function to save a figure
function saveFigure(filename, figNumber)
  exportgraphics(figure(figNumber), filename);
end

function A = psf2A(psf, bc)
  assert(size(psf, 1) == size(psf, 2));
  m = size(psf, 1);
  nzmax = m * m * nnz(psf);
  A_i = zeros(nzmax, 1);
  A_j = zeros(nzmax, 1);
  A_s = zeros(nzmax, 1);
  idx = 1;
  for i = 1:m * m
      x = zeros(m, m);
      x(i) = 1;
      x_conv_psf = imfilter(x, psf, 'conv', 'same', bc);
      [row, col, val] = find(x_conv_psf);
      numVals = length(val);
      A_i(idx:idx+numVals-1) = sub2ind([m, m], row, col);
      A_j(idx:idx+numVals-1) = i;
      A_s(idx:idx+numVals-1) = val;
      idx = idx + numVals;
  end
  A_i = A_i(1:idx-1);
  A_j = A_j(1:idx-1);
  A_s = A_s(1:idx-1);
  A = sparse(A_i, A_j, A_s, m * m, m * m);
end