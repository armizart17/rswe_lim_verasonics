function [fitresult, gof] =createFit_atenuation_1d_RSWE (x_v, R_cc_1d)
%CREATEFIT(X_V,R_CC_1D)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : x_v
%      Y Output: R_cc_1d
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 09-Jan-2023 19:42:55


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( x_v, R_cc_1d );

% Set up fittype and options.
ft = fittype( '3*((-k)*x*cos(k*x) + sin(k*x) + k^2*x^2*sin(k*x))/(4*k^3*x^3)', 'independent', 'x', 'dependent', 'y' );
excludedPoints = yData >= 1;
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = 50;
opts.StartPoint = 0.981728773776423;
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
% 
% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData, excludedPoints );
% legend( h, 'R_cc_1d vs. x_v', 'Excluded R_cc_1d vs. x_v', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'x_v', 'Interpreter', 'none' );
% ylabel( 'R_cc_1d', 'Interpreter', 'none' );
% grid on


