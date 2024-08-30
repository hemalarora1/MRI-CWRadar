% From MSSA_beginners_guide_v3, modify initial waveforms to test mSSA
clear;    % clear workspace

% Set general Parameters
M = 30;    % window length of SSA
N = 200;   % length of generated time series
T = 22;    % period length of sine function
stdnoise = 0.1; % noise-to-signal ratio

% Create time series X
% First of all, we generate two time series, a sine function
% of length N and the same function to the power of 3, both
% with observational white noise.

t = (1:N)';
X1 = sin(2*pi*t/T);     % sine function
X2 = cos(2*pi*t/T).^3;  % nonlinear transformation
X3 = sin(2*pi*3.5*t/T);
X4 = cos(2*pi*1.8*t/T).^3;
%X3 = sin(2*pi*t/T).^5;     % sine function
%X4 = cos(2*pi*t/T).^7;
noise = stdnoise*randn(N,4);  % Gaussian noise
X1 = X1+noise(:,1);
X2 = X2+noise(:,2);
X3 = X3+noise(:,3);
X4 = X4+noise(:,4);
X1 = X1-mean(X1); % remove mean value
X2 = X2-mean(X2);
X3 = X3-mean(X3);
X4 = X4-mean(X4);
X1 = X1/std(X1);  % normalize to std=1
X2 = X2/std(X2);
X3 = X3/std(X3);
X4 = X4/std(X4);

X = X1;
%X = [X1 X2]; % multivariate time series
%X = [X1 X2 X3];
X = [X1 X2 X3 X4];

% if M=1, should degenerate to PCA, but causes error in RC indices
% as there is no diagonal term? Also if M=1, by dfn only 1 RC so
% range of RC terms is only single term.

[RC, LAMBDA] = mSSA(X, M, 1);
Xrec = recMSSA(RC, 1, 2, 1) ;
