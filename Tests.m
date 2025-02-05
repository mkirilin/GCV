% TODO: 1. find optimal SNRs
%         (a) find optimal SNR for blurGauss
%           (i)   128 -- (-2:2)
%           (ii)  256 -- (-2:2)
%           (iii) 512 -- (-2:2)
%         (b) find optimal SNR for CT
%           (i)   128 -- (-2:2)
%           (ii)  256 -- (-2:2)
%**************************************
% Clear workspace and window
clear; clc;

% Choose if you would like to see the results displayed in a single figure 
% window ('subplots') or in multiple figure windows ('manyplots')
dispres = 'manyplots';

LW = 2;  % Plot line width
MS = 10; % Size of markers on plots

numTests = 10;
SNR = 10.^(-2:2)';

model = 'blurGauss';  % Choose between 'blur', 'blurGauss' and 'CT'
n_values = [64];  % List of n values to iterate over

resultsDir = fullfile(fileparts(mfilename('fullpath')), 'Results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

for n = n_values
tic
  % Define the test problem
  [A, b, x, ProbInfo, m0] = defineTestProblem(model, n);

  m = min(size(A)); % Maximum possible rank of A
  if strcmp(model, 'CT')
    m_sv = sprank(A); % Structural rank of A (>= rank(A))
  else
    m_sv = m;
  end
  %assert(m0 <= m_sv, 'm0 must be less than or equal to the number of rows of A');
  %assert (sprank(A) == m_sv, 'A must have full rank');

  if strcmp(model, 'blurGauss')
    [~, S, ~] = svds(A, b, m_sv);

    % Sort singular values in descending order and change U and V accordingly
    [S, idx] = sort(S, 'descend');

  else
    [U,S,V] = svds(A, m_sv);
    S = diag(S);
  end

  % Plot singular values decay of A
  figure(5); clf;
  semilogy(S, 'linewidth', LW);
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
      if strcmp(model, 'blurGauss')
        [X_gcv, X_opt, error_gcv, error_opt, k_gcv, k_opt] =...
          gcvBlurGauss(S(1:m0), x, bn, m0, m, false, idx, n);
      else
        [X_gcv, X_opt, error_gcv, error_opt, k_gcv, k_opt] =...
          gcv(U(:,1:m_sv), S(1:m_sv), V(:,1:m_sv), x, bn, m_sv, m, false);
      end
      testData = [testData; table(i, error_gcv, "gcv",...
                                  'VariableNames', {'SNR', 'Error', 'Method'})];
      testData = [testData; table(i, error_opt, "opt",...
                                  'VariableNames', {'SNR', 'Error', 'Method'})];
      kData = [kData; table(i, k_gcv, "gcv",...
                            'VariableNames', {'SNR', 'k', 'Method'})];
      kData = [kData; table(i, k_opt, "opt",...
                            'VariableNames', {'SNR', 'k', 'Method'})];
    end
  end
  toc

  % Plot results
  plotErrorResults(testData, n, m0, m, SNR);
  plotKResults(kData, n, m0, m, SNR);

  % Display the reconstructions
  displayReconstructions(x, bn, X_gcv, X_opt, ProbInfo);

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
    options.CommitCrime = 'on';
    m0 = n*n;
    [A, b, x, ProbInfo] = PRtomo(n, options);
  elseif strcmp(model, 'blurGauss')
    options = IRset();
    %options.BlurLevel = 'severe';
    %options.BC = 'zero';
    %options.trueImage = 'satellite';
    %options.CommitCrime = 'on';
    m0 = n*n;
    [A, b, x, ProbInfo] = PRblurgauss(n, options);
  else
    error('Invalid model');
  end
end

% Function to plot results
function plotErrorResults(testData, n, m0, m, SNR)
  figure(10); clf;
  title(['Results for n = ', num2str(n), ', ', num2str(m0), '/', num2str(m),...
        ' singular values'], 'interpreter', 'latex', 'fontsize', 18);
  boxchart(testData.SNR, testData.Error, 'GroupByColor', testData.Method);
  
  xlabel('SNR');
  xtick=SNR;
  xticklab = cellstr(num2str(round(log10(xtick(:))), '10^{%d}'));
  set(gca,'XTickLabel',xticklab,'TickLabelInterpreter','tex')
  ax = gca;
  ylabel('Relative error');
  ax.YAxis.Scale ="log";
  legend()
end

function plotKResults(kData, n, m0, m, SNR)
  figure(11); clf;
  title(['Results for n = ', num2str(n), ', ', num2str(m0), '/', num2str(m), ' singular values'],...
        'interpreter', 'latex', 'fontsize', 18);
  boxchart(kData.SNR, kData.k, 'GroupByColor', kData.Method);

  xlabel('SNR');
  xtick=SNR;
  xticklab = cellstr(num2str(round(log10(xtick(:))), '10^{%d}'));
  set(gca,'XTickLabel',xticklab,'TickLabelInterpreter','tex')
  %ax = gca;
  ylabel('k');
  %ax.YAxis.Scale ="log";
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
  savefig(figure(figNumber), [filename(1:end-3), 'fig']);
end
