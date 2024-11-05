
function [RC,LAMBDA] = compSSA(X, WinSiz, ShowFigs, saveDir, position)
%compSSA performs singular spectral analysis on complex time series Z. 
%
%   X: complex time series (has both real and imaginary
%   WinSiz : integer represent fraction of length to use as window
%   ShowFigs : 0 don't show figs, 1 show figs
%
% example:
%
%   [RC,LAMBDA] = compSSA(Z, 4, 1) ;
%  causes the SSA sliding window to be 1/4 of the input length N
%  so WinSiz=4 means the sliding window M=round(N/4). Figures will be
%  shown.
%  


%Z = I+j*Q;
%Z = Z-mean(Z);  % subtracting mean 
N = length(X);
M = round(N/WinSiz);

% Choose covariance estimation
%Choose between Toeplitz approach (cf. Vautard & Ghil) and trajectory
%approach (cf. Broomhead & King).

if(0)
% Calculate covariance matrix C (Toeplitz approach)
%Next, we calculate the covariance matrix. There are several numerical
%approaches to estimate C. Here, we calculate the covariance function with
%CORR and build C with the function TOEPLITZ.
covX = xcorr(X,M-1,'unbiased');
Ctoep=toeplitz(covX(M:end));
C = Ctoep;
end

if (1)
% Calculate covariance matrix (trajectory approach)
%An alternative approach is to determine C directly from the scalar product
%of Y, the time-delayed embedding of X. Although this estimation of C does
%not give a Toeplitz structure, with the eigenvectors not being symmetric 
%or antisymmetric, it ensures a positive semi-definite covariance matrix.
Y=zeros(N-M+1,M);
for m=1:M
  Y(:,m) = X((1:N-M+1)+m-1);
end;
Cemb=Y'*Y / (N-M+1);
C=Cemb;
%size(C)
end

% Calculate eigenvalues LAMBDA and eigenvectors RHO
%In order to determine the eigenvalues and eigenvectors of C, we use the
%function EIG. This function returns two matrices, the matrix RHO with 
%eigenvectors arranged in columns, and the matrix LAMBDA with eigenvalues
%along the diagonal.
[RHO,LAMBDA] = eig(C);
LAMBDA = diag(LAMBDA);               % extract the diagonal elements
[LAMBDA,ind]=sort(LAMBDA,'descend'); % sort eigenvalues
RHO = RHO(:,ind);                    % and eigenvectors

% Calculate principal components PC
%The principal components are given as the scalar product between Y,
%the time-delayed embedding of X, and the eigenvectors RHO.
PC = Y*RHO; % principal components

% Calculate reconstructed components RC
%In order to determine the reconstructed components RC, we have to invert
%the projecting PC = Y*RHO; i.e. RC = Y*RHO*RHO'=PC*RHO'. Averaging along
%anti-diagonals gives the RCs for the original input X.
RC=zeros(N,M); % reconstruct components RC
for m=1:M
  buf=PC(:,m)*RHO(:,m)'; % invert projection
  buf=buf(end:-1:1,:);
  for n=1:N % anti-diagonal averaging
    RC(n,m)=mean( diag(buf,-(N-M+1)+n) );
  end
end;

if(ShowFigs)
t = 1:length(X);
figure(1);
set(gcf,'name','Time series X');
clf;
plot(t,real(X),'b-');
hold on;
plot(t,imag(X),'r-');
hold off;
saveas(gcf, fullfile(saveDir, sprintf('timeseries_X_%s.png', position)));

figure(2);
set(gcf,'name','Covariance matrix');
clf;
imagesc(abs(C));
%imagesc(Cemb);
axis square
%set(gca,'clim',[-1 1]);
colorbar
saveas(gcf, fullfile(saveDir, sprintf('covariance_matrix_%s.png', position)));

figure(3);
set(gcf,'name','Eigenvalues LAMBDA')
clf;
semilogy(LAMBDA,'o-');
saveas(gcf, fullfile(saveDir, sprintf('eigenvalues_lambda_%s.png', position)));

  if(0)
    figure(4);
    set(gcf,'name','Principal components PCs')
    clf;
    for m=1:4
        subplot(4,1,m);
        plot(t(1:N-M+1),PC(:,m),'k-');
        ylabel(sprintf('PC %d',m));
        % ylim([-10 10]);
    end;
    saveas(gcf, fullfile(saveDir, sprintf('principal_components_PCs_%s.png', position)));
  end
  
if(1)
figure(5);
set(gcf,'name','Reconstructed components RCs')
clf;
for m=1:6
  subplot(6,1,m);
  %plot(t,RC(:,m),'r-');
  plot(t,real(RC(:,m)),'b-',t,imag(RC(:,m)),'r-');
  ylabel(sprintf('RC %d',m));
%  ylim([-1 1]);
end;
saveas(gcf, fullfile(saveDir, sprintf('Reconstructed_components_RCs_%s.png', position)));
end

% Compare reconstruction and original time series
%Note that the original time series X can be completely reconstructed by
%the sum of all reconstructed components RC (upper panel). The sine
%function can be reconstructed with the first pair of RCs (lower panel).
figure(6);
NS=1; NP = 30;
set(gcf,'name','Original time series X and reconstruction RC')
clf;
subplot(2,1,1);
%plot(t,X,'b-',t,sum(RC(:,:),2),'r-');
%legend('Original','Complete reconstruction');
plot(t,real(X),'b-',t,real(sum(RC(:,NS:NP),2)),'r-');
legend('Original','Reconstruction with RCs');

subplot(2,1,2);
plot(t,imag(X),'b-',t,imag(sum(RC(:,NS:NP),2)),'r-');
legend('Original','Reconstruction with RCs');
saveas(gcf, fullfile(saveDir, sprintf('timeseriesX_and_reconstructionRCs_%s.png', position)));

if(0)
figure(7);
NS=1; NP=30;
rI = real(sum(RC(:,NS:NP),2));
rQ = imag(sum(RC(:,NS:NP),2));

%plot(t,real(X),'b-',t,real(sum(RC(:,NS:NP),2)),'r-');
plot(t,real(sum(RC(:,NS:NP),2)),'b-'); % sum NS:NP RCS
hold on;
%plot(t,imag(X),'k-',t,imag(sum(RC(:,NS:NP),2)),'m-');
plot(t,imag(sum(RC(:,NS:NP),2)),'r-'); % sum NS:NP RCs
hold off;
saveas(gcf, fullfile(saveDir, sprintf('unknownfig7_compSSA%s.png', position)));
end
end  % ShowFigs

end



