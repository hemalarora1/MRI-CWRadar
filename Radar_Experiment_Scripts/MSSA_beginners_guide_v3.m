% M-SSA tutorial with Matlab
% This Matlab tutorial demonstrates step by step the multivariate
% singular spectrum analysis. The steps are almost similar
% to those of a singular spectrum analysis.

% Copyright (c) 2013-2016, Andreas Groth, University of California, Los Angeles.
% All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% - Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% - Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

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
%X3 = sin(2*pi*3.5*t/T);
%X4 = cos(2*pi*1.8*t/T).^3;
X3 = sin(2*pi*t/T).^5;     % sine function
X4 = cos(2*pi*t/T).^7;
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

X = [X1 X2]; % multivariate time series

figure(1);
clf;
set(1,'name','Time series X');
subplot(1,2,1);
plot(X1, 'r-');
title('Time series X1');
subplot(1,2,2);
plot(X2, 'r-');
title('Time series X2');

% Calculate covariance matrix C (Toeplitz approach)
% Next, we calculate the covariance matrix. There are several numerical
% approaches to estimate C. Here, we calculate the covariance function with
% CORR and build C with the function TOEPLITZ.

if(0)
covXX=xcorr(X1,M-1,'unbiased');
covYY=xcorr(X2,M-1,'unbiased');
covXY = xcorr(X1,X2,M-1,'unbiased');

C11=toeplitz(covXX(M:end));
C21=toeplitz(covXY(M:-1:1),covXY(M:end));
C12=C21';
C22=toeplitz(covYY(M:end));

Ctoep = [C11 C12 ;...
         C21 C22  ];
       
figure(2);
set(gcf,'name','Covariance matrix');
clf;
imagesc(Ctoep);
axis square
set(gca,'clim',[-1 1]);
colorbar
end


% Calculate covariance matrix (trajectory approach)
% This alternative approach is to determine C directly from the scalar
% product of Y, the time-delayed embedding of X. Although this estimation
% of C does not give a Toeplitz structure, with the eigenvectors not being
% symmetric or antisymmetric, it ensures a positive semi-definite covariance
% matrix.

Y1=zeros(N-M+1,M);
Y2=zeros(N-M+1,M);
for m=1:M                 % create time-delayed embedding of X
  Y1(:,m) = X1((1:N-M+1)+m-1);
  Y2(:,m) = X2((1:N-M+1)+m-1);
end;
Y = [Y1 Y2];
Cemb=Y'*Y / (N-M+1);

figure(2);
set(gcf,'name','Trajectory Matrices');
clf;
imagesc(Y);
set(gca,'clim',[-1 1]);
colorbar

figure(3);
set(gcf,'name','Covariance matrix');
clf;
imagesc(Cemb);
axis square
set(gca,'clim',[-1 1]);
colorbar

%% Choose covariance estimation
% Choose between Toeplitz approach (cf. Vautard & Ghil) and trajectory
% approach (cf. Broomhead & King).

% C=Ctoep;
C=Cemb;


%% Calculate eigenvalues LAMBDA and eigenvectors RHO
% In order to determine the eigenvalues and eigenvectors of C
% we use the function EIG. This function returns two matrices,
% the matrix RHO with eigenvectors arranged in columns, and the
% matrix LAMBDA with eigenvalues on the diagonal.

[RHO,LAMBDA] = eig(C);
LAMBDA = diag(LAMBDA);               % extract the diagonal
[LAMBDA,ind]=sort(LAMBDA,'descend'); % sort eigenvalues
RHO = RHO(:,ind);                    % and eigenvectors

figure(4);
clf;
set(gcf,'name','Eigenvectors RHO and eigenvalues LAMBDA')
subplot(3,1,1);
plot(LAMBDA,'o-');
subplot(3,1,2);
plot(RHO(:,1:2), '-');
legend('1', '2');
subplot(3,1,3);
plot(RHO(:,3:4), '-');
legend('3', '4');

%% Calculate principal components PC
% The principal components are given as the scalar product
% between Y, the time-delayed embedding of X1 and X2, and the
% eigenvectors RHO.

PC = Y*RHO;

figure(5);
set(gcf,'name','Principal components PCs')
clf;
for m=1:4
  subplot(4,1,m);
  plot(PC(:,m),'k-');
  ylabel(sprintf('PC %d',m));
  ylim([-10 10]);
end;

%% Calculate reconstructed components RC
% In order to determine the reconstructed components RC, we have to invert
% the projecting PC = Y*RHO; i.e. RC = Y*RHO*RHO'=PC*RHO'. Averaging along
% anti-diagonals gives the RCs for the original input X.

RC1=zeros(N,2*M);
RC2=zeros(N,2*M);
for m=1:2*M
  buf1=PC(:,m)*RHO(1:M,m)'; % invert projection - first channel
  buf1=buf1(end:-1:1,:);
  
  buf2=PC(:,m)*RHO(M+1:end,m)'; % invert projection - second channel
  buf2=buf2(end:-1:1,:);
  
  for n=1:N % anti-diagonal averaging
    RC1(n,m)=mean( diag(buf1,-(N-M+1)+n) );  
    RC2(n,m)=mean( diag(buf2,-(N-M+1)+n) );  
  end
end;

figure(6);
set(gcf,'name','Reconstructed components RCs')
clf;
for m=1:4
  subplot(4,2,2*m-1);
  plot(RC1(:,m),'r-');
  ylabel(sprintf('RC %d',m));
  ylim([-1 1]);
  
  subplot(4,2,2*m);
  plot(RC2(:,m),'r-');
  ylabel(sprintf('RC %d',m));
  ylim([-1 1]);
end;

%% Compare reconstruction and original time series
% Note that the original time series X can be completely reconstructed by
% the sum of all reconstructed components RC (upper panels). The sine
% function (lower left panel) can be reconstructed with the first pair of
% RCs, where more components are required for the nonlinear oscillation
% (lower right panel).

figure(7);
set(gcf,'name','Original time series X and reconstruction RC')
clf;
subplot(2,2,1)
plot(t,X1,'b-',t,sum(RC1,2),'r-');
subplot(2,2,2)
plot(t,X2,'b-',t,sum(RC2,2),'r-');
legend('Original','full reconstruction');

subplot(2,2,3)
plot(t,X1,'b',t,sum(RC1(:,1:2),2),'r');
subplot(2,2,4)
plot(t,X2,'b',t,sum(RC2(:,1:2),2),'r');
legend('Original','RCs 1-2');
