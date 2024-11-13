% EXblur_cgls_hybrid Example script, speckle deblurring problem

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% clear workspace and window
clear, clc

% Choose if you would like to see the results displayed in a single figure 
% window ('subplots') or in multiple figure windows ('manyplots')
% dispres = 'subplots';
dispres = 'manyplots';

LW = 2;  % Plot line width
MS = 10; % Size of markers on plots

rng(0);  % Make sure this test is repeatable.

% Define the test problem.
NoiseLevel = 0.01;
n = 64;
[A, b, x, ProbInfo] = PRblurrotation(n);
[bn, NoiseInfo] = PRnoise(b, 'gauss', NoiseLevel);
[X_gcv, X_opt] = gcv(A, x, bn, 100);

% Display the reconstructions;
% uncomment as appropriate to avoid displaying titles and legends
strcmp(dispres, 'manyplots')
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
    %figure(3), clf
    %axes('FontSize', 24), hold on
    %semilogy(1:100, IterInfo_cgls.Enrm, 'b-', 'LineWidth', LW)
    %hold on
    %semilogy(0:100, [norm(bn); IterInfo_cgls.Enrm], 'k-.', 'LineWidth', LW)
    %semilogy(IterInfo_cgls.BestReg.It, IterInfo_cgls.BestReg.Enrm, 'ro', 'LineWidth', LW, 'MarkerSize', MS)
    %semilogy(IterInfo_cgls_dp.its, IterInfo_cgls_dp.Enrm(end), 'ms', 'LineWidth', LW, 'MarkerSize', MS)
    %hl = legend('{\tt IRcgls} errors','{\tt IRhybrid\_lsqr} errors', ...
    %  'optimal {\tt IRcgls} stopping iteration','{\tt IRcgls} DP stopping iteration', ...
    %  '{\tt IRhybrid\_lsqr} DP stopping iteration');
    %set(hl,'interpreter','latex','fontsize',18)
    % title('Error history','interpreter','latex','fontsize',18)
    %axis([0,100,0.15,IterInfo_cgls.Enrm(1)])
    %
    figure(3), clf
    PRshowx(X_gcv, ProbInfo)
    title('GCV sol.',...
    'interpreter','latex','fontsize',18)
    figure(4), clf
    PRshowx(X_opt, ProbInfo)
    title('Optimal sol.',...
    'interpreter','latex','fontsize',18)

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
    exportgraphics(figure(1), 'EXblur_cgls_hybrid_a.eps', 'ContentType', ...
                 'vector', 'Resolution', 300);
    exportgraphics(figure(2), 'EXblur_cgls_hybrid_b.eps', 'ContentType', ...
                 'vector', 'Resolution', 300);
    exportgraphics(figure(3), 'EXblur_cgls_hybrid_c.eps', 'ContentType', ...
                 'vector', 'Resolution', 300);
    exportgraphics(figure(4), 'EXblur_cgls_hybrid_d.eps', 'ContentType', ...
                 'vector', 'Resolution', 300);
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
    saveas(figure(1), 'EXblur_cgls_hybrid_a.fig')
    saveas(figure(2), 'EXblur_cgls_hybrid_b.fig')
    saveas(figure(3), 'EXblur_cgls_hybrid_c.fig')
    saveas(figure(4), 'EXblur_cgls_hybrid_d.fig')
end
cd(oldcd)

n = 64;
NoiseLevel = 0.1;
options = IRset();
options.sm = true;
[A, b, x, ProbInfo] = PRtomo(n, options);
[bn, NoiseInfo] = PRnoise(b, 'gauss', NoiseLevel);
[X_gcv_tomo, X_opt_tomo] = gcv(A, x, bn, 4096);


figure(10), clf
PRshowx(x, ProbInfo)
set(gca,'fontsize',24)
title('True solution','interpreter','latex','fontsize',18)

figure(20), clf
PRshowb(b, ProbInfo)
set(gca,'fontsize',24)
title('Noisy data','interpreter','latex','fontsize',18)

figure(30), clf
PRshowx(X_gcv_tomo, ProbInfo)
set(gca,'fontsize',24)
title('GCV sol.','interpreter','latex','fontsize',18)

figure(40), clf
PRshowx(X_opt_tomo, ProbInfo)
set(gca,'fontsize',24)
title('Opt sol.','interpreter','latex','fontsize',18)