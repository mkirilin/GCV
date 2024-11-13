% clear workspace and window
clear, clc

% Choose if you would like to see the results displayed in a single figure 
% window ('subplots') or in multiple figure windows ('manyplots')
% dispres = 'subplots';
dispres = 'manyplots';

LW = 2;  % Plot line width
MS = 10; % Size of markers on plots

rng(0);  % Make sure this test is repeatable.
model = 'CT';
NoiseLevel = 0.1;
n = 64;

if strcmp(model, 'blur')
  % Define the test problem.
  [A, b, x, ProbInfo] = PRblurrotation(n);
  [bn, NoiseInfo] = PRnoise(b, 'gauss', NoiseLevel);
  [X_gcv, X_opt] = gcv(A, x, bn, 100);
elseif strcmp(model, 'CT')
  n = 64;
  options = IRset();
  options.sm = true;
  [A, b, x, ProbInfo] = PRtomo(n, options);
  [bn, NoiseInfo] = PRnoise(b, 'gauss', NoiseLevel);
  [X_gcv, X_opt] = gcv(A, x, bn, 4096);
end

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
