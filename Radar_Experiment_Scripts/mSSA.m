% this function adapted from code MSSA_beginners_guide_v3.m written by Andreas GrothCopyright (c) 2013-2016, Andreas Groth, University of California, Los Angeles.
% All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% - Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% - Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

function [RC,LAMBDA] = mSSA(X, M, ShowFigs)
%mSSA performs multi-channel singular spectral analysis on complex time
%     array X of several time series.  This code will determine the
%     number of time series from the second dimension size of X. The
%     length of a time series will be determined from the first dimension
%
%   X: array of complex time series (has both real and imaginary
%   M : integer represent  length to use as window trajectory
%   ShowFigs : 0 don't show figs, 1 show figs
%
% example:
%
%   [RC,LAMBDA] = mSSA([X1, X2, X3], 4, 1) ;
%  causes the SSA sliding window to be 1/4 of the input length N of
%  X1, X2, or X3 so WinSiz=4 means the sliding window M=round(N/4).
%  Figures will be shown.  Here, 3 "channels" are input into mSSA
%  For now, assume input array may have been prewhitened by user
%  ie zero mean, standard deviation of 1.
%
%Z = I+j*Q;
%Z = Z-mean(Z);  % subtracting mean 


dimX = size(X);
N = dimX(1) ;   % length of each time series - must be equal, not checked
nCH = dimX(2);  % number of channels

% Covariance estimation:
% Calculate covariance matrix (trajectory approach) of Broomhead & King.
%This alternative approach is to determine C directly from the scalar product
%of Y, the time-delayed embedding of X. Although this estimation of C does
%not give a Toeplitz structure, with the eigenvectors not being symmetric 
%or antisymmetric, it ensures a positive semi-definite covariance matrix.

if(1)  % from varimax tutorial, created trajectory
    D=nCH;
    xtde=zeros(N,N,D);
    for d=1:D
        xtde(:,:,d)=hankel(X(:,d));
    end
    xtde=xtde(1:N-M+1,1:M,:);
    xtde=reshape(xtde,N-M+1,D*M,1);
    Y = xtde;
    % M-SSA analysis
    %C=xtde'*xtde/(N-M+1); % Broomhead and King (1986)
end

if(0)
    Yc = zeros(N-M+1,M);
    Y = [];  % null matrix will grow
    for ch = 1:nCH
        Xc = X(:,ch) ;
        for m=1:M                 % create time-delayed embedding of X
            Yc(:,m) = Xc((1:N-M+1)+m-1);
        end;
        Y = [Y Yc];
    end
end

C=Y'*Y / (N-M+1);  % Broomhead and King (1986)
Ysiz = size(Y)
Csiz = size(C)

% Calculate eigenvalues LAMBDA and eigenvectors RHO
%In order to determine the eigenvalues and eigenvectors of C, we use the
%function EIG. This function returns two matrices, the matrix RHO with 
%eigenvectors arranged in columns, and the matrix LAMBDA with eigenvalues
%along the diagonal.
[RHO,LAMBDA] = eig(C);
LAMBDA = diag(LAMBDA);               % extract the diagonal elements
[LAMBDA,ind]=sort(LAMBDA,'descend'); % sort eigenvalues
RHO = RHO(:,ind);                    % and eigenvectors
RHOsiz = size(RHO)
% Calculate principal components PC
%The principal components are given as the scalar product between Y,
%the time-delayed embedding of X1, X2 .. Xnch, and the eigenvectors RHO.
PC = Y*RHO; % principal components
%figure(1);
%imagesc(PC);
PCsiz = size(PC)

% Calculate reconstructed components RC
% In order to determine the reconstructed components RC, we have to invert
% the projecting PC = Y*RHO; i.e. RC = Y*RHO*RHO'=PC*RHO'. Averaging along
% anti-diagonals gives the RCs for the original input X.

% why RC now nCH*M wide? really
RC = zeros(N, nCH*M, nCH); % reconstruct components RC.  
for m = 1:nCH*M
  for nc = 1:nCH
    Mo = (nc-1)*M+1; Ms = nc*M;  % get range start/stop index
    buf=PC(:,m)*RHO(Mo:Ms,m)'; % invert projection of channel CH
    buf=buf(end:-1:1,:);
    for n=1:N % anti-diagonal averaging
        RC(n,m, nc) = mean( diag(buf, -(N-M+1)+n) );  
    end
  end
end

if(ShowFigs)
     % ylabel(sprintf('PC %d',m));
t = 1:N;
if(1)
figure(11);
set(gcf,'name','Time series X');
clf;
for m=1:nCH
   subplot(1,nCH,m);
   plot(real(X(:,m)), 'b-');
   hold on;
   plot(imag(X(:,m)), 'r-');
   hold off;
   label = sprintf('Time series X%d',m);
   title(label);
end
end

if(1)
figure(12);
set(gcf,'name','Trajectory Matrices');
clf;
imagesc(Y);
set(gca,'clim',[-1 1]);
colorbar

figure(13);
set(gcf,'name','Covariance matrix');
clf;
%imagesc(abs(C));
%imagesc(real(C));
imagesc(C);
axis square
set(gca,'clim',[-1 1]);
colorbar
end

if(1)
figure(14);
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
end

if(1)
figure(15);
set(gcf,'name','Principal components PCs')
clf;
for m=1:4
  subplot(4,1,m);
  plot(PC(:,m),'k-');
  ylabel(sprintf('PC %d',m));
  ylim([-10 10]);
end;
end

if(1)
figure(16);
set(gcf,'name','Reconstructed components RCs')
clf;
for m=1:4
 for k=1:nCH
  subplot(4,nCH,nCH*(m-1)+k);
  plot(t,real(RC(:,m,k)),'b-',t,imag(RC(:,m,k)),'r-');
  ylabel(sprintf('RC %d',m));
  ylim([-1 1]);
 end
end
end


% Compare reconstruction and original time series
%Note that the original time series X can be completely reconstructed by
%the sum of all reconstructed components RC (upper panel). The sine
%function can be reconstructed with the first pair of RCs (lower panel).
if(1)
figure(17);
set(gcf,'name','Original time series X and reconstruction RC')
clf;
for m=1:nCH
    subplot(2,nCH,m)
    plot(t,real(X(:,m)),'b-',t,real(sum(RC(:,:,m),2)),'r-');
end
legend('Original','full reconstruction');
for m=1:nCH
    subplot(2,nCH,m+nCH)
    plot(t,real(X(:,m)),'b',t,real(sum(RC(:,1:2,m),2)),'r');   
end
legend('Original','RCs 1-2');
end

% next lines old code. Not fixed for mSSA yet.
if(0)
figure(7);
NS=1; NP=10;
rI = real(sum(RC(:,NS:NP),2));
rQ = imag(sum(RC(:,NS:NP),2));

%plot(t,real(X),'b-',t,real(sum(RC(:,NS:NP),2)),'r-');
plot(t,real(sum(RC(:,NS:NP),2)),'b-'); % sum NS:NP RCS
hold on;
%plot(t,imag(X),'k-',t,imag(sum(RC(:,NS:NP),2)),'m-');
plot(t,imag(sum(RC(:,NS:NP),2)),'r-'); % sum NS:NP RCs
hold off;
end
end  % ShowFigs

end



