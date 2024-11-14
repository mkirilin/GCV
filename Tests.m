% clear workspace and window
clear, clc

% Choose if you would like to see the results displayed in a single figure 
% window ('subplots') or in multiple figure windows ('manyplots')
% dispres = 'subplots';
dispres = 'manyplots';

LW = 2;  % Plot line width
MS = 10; % Size of markers on plots

numTests = 5;
SNR = 10.^(0:4)';

model = 'CT';  % Choose between 'blur' and 'CT'
n = 64;

% Preallocate table for test data
testData = table('Size', [0, 3], 'VariableTypes', {'double', 'double', 'string'},...
                 'VariableNames', {'SNR', 'Error', 'Method'});

% Define the test problem.
if strcmp(model, 'blur')
  m0 = 100;
  [A, b, x, ProbInfo] = PRblurrotation(n);
elseif strcmp(model, 'CT')
  options = IRset();
  options.sm = true;
  m0 = 200;
  [A, b, x, ProbInfo] = PRtomo(n, options);
end

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

figure(10), clf;
testData.SNR = categorical(testData.SNR, SNR, {'1', '1e1', '1e2', '1e3', '1e4'});
boxchart(testData.SNR, testData.Error, 'GroupByColor', testData.Method);

return
% Display the reconstructions;
% uncomment as appropriate to avoid displaying titles and legends
figure(1), clf
PRshowx(x, ProbInfo)
set(gca,'fontsize',24)
title('True solution','interpreter','latex','fontsize',18)
%
figure(2), clf
PRshowb(b, ProbInfo)
set(gca,'fontsize',24)
title('Noisy data','interpreter','latex','fontsize',18)
%
figure(3), clf
PRshowx(X_gcv, ProbInfo)
title('GCV sol.','interpreter','latex','fontsize',18)
%
figure(4), clf
PRshowx(X_opt, ProbInfo)
title('Optimal sol.','interpreter','latex','fontsize',18)

%return

% A number of instructions useful to save the displayed figures follow;
% the default is not to execute them. If you wish to save the displayed
% figures in the dedicated 'Results' folder, please comment the above
% return statement

currentFolder = fileparts(mfilename('fullpath'));
cd(currentFolder)
oldcd = cd;
if strcmp(dispres, 'subplots')
  try
    cd('Results')
  catch
    mkdir('Results')
    cd('Results')
  end
  figure(1), print -dpng -r300 EXblur_cgls_hybrid
elseif strcmp(dispres, 'manyplots')
  try
    cd('Results')
  catch
    mkdir('Results')
    cd('Results')
  end
  exportgraphics(figure(1), 'Orig.eps', 'ContentType', 'vector', 'Resolution', 300);
  exportgraphics(figure(2), 'Signal.eps', 'ContentType', 'vector', 'Resolution', 300);
  exportgraphics(figure(3), 'GCV.eps', 'ContentType', 'vector', 'Resolution', 300);
  exportgraphics(figure(4), 'Opt.eps', 'ContentType', 'vector', 'Resolution', 300);

  %saveas(figure(1), 'Orig.fig', 'fig');
  %saveas(figure(2), 'Signal.fig');
  %saveas(figure(3), 'GCV.fig');
  %saveas(figure(4), 'Opt.fig');
end
cd(oldcd)
