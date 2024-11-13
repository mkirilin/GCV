% clear workspace and window
clear, clc

% Choose if you would like to see the results displayed in a single figure 
% window ('subplots') or in multiple figure windows ('manyplots')
% dispres = 'subplots';
dispres = 'manyplots';

LW = 2;  % Plot line width
MS = 10; % Size of markers on plots

numTests = 2;
errors_gcv = zeros(numTests,1);
errors_opt = zeros(numTests,1);

model = 'CT';
NoiseLevel = 0.1;
n = 64;

for i = 1:numTests
  fprintf('Test %d\n', i);
  rng(i);  % Set seed for reproducibility

  if strcmp(model, 'blur')
    % Define the test problem.
    [A, b, x, ProbInfo] = PRblurrotation(n);
    [bn, NoiseInfo] = PRnoise(b, 'gauss', NoiseLevel);
    [X_gcv, X_opt] = gcv(A, x, bn, 100);
  elseif strcmp(model, 'CT')
    options = IRset();
    options.sm = true;
    [A, b, x, ProbInfo] = PRtomo(n, options);
    NoiseLevel = norm(b) / sqrt(size(b,1));
    [bn, NoiseInfo] = PRnoise(b, 'gauss', NoiseLevel);
    [X_gcv, X_opt] = gcv(A, x, bn, 4096);
  end
  errors_gcv(i) = norm(X_gcv - x);
  errors_opt(i) = norm(X_opt - x);
end

%figure(10), clf;
%boxplot([errors_gcv, errors_opt], {'GCV','Optimal'});
%title('Error Distribution over 10 Runs');
% Plot a boxplot of errors_gcv in blue and errors_opt in red
figure(10), clf;
boxplot(errors_gcv, 'Color', 'b');
hold on;
boxplot(errors_opt, 'Color', 'r');
title('Error Distribution over 10 Runs');
xlabel('Method');
ylabel('Error');
legend('GCV', 'Optimal');

%return
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
end
cd(oldcd)

% Uncomment the following return statement if you wish to save the
% displayed figures as MATLAB figures

% return

oldcd = cd;
if strcmp(dispres, 'subplots')
  try
    cd('Results')
  catch
    mkdir('Results')
    cd('Results')
  end
  figure(1), saveas('EXblur_cgls_hybrid.fig')
elseif strcmp(dispres, 'manyplots')
  try
    cd('Results')
  catch
    mkdir('Results')
    cd('Results')
  end
  saveas(figure(1), 'Orig.fig')
  saveas(figure(2), 'Signal.fig')
  saveas(figure(3), 'GCV.fig')
  saveas(figure(4), 'Opt.fig')
end
cd(oldcd)
