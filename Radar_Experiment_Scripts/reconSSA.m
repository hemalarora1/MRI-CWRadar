function Zssa = reconSSA(RC, NS, NP, ShowFigs, saveDir, position )
% reconstruct time series from SSA components NS through NP.
% example:
%    Drecon = reconSSA(RC, 2, 10, 0) ;
%    uses 2nd through 10th SSA componnents. No figure shown ,
%
%   Drecon = reconSSA(RC, 1, 20, 1) ;
%    uses first 20 SSA components to reconstruct data. Figures are shown.
%

Zssa = (sum(RC(:,NS:NP),2));
Xo = sum(RC(:,:),2);  % sum all terms to recover original data

if(ShowFigs)
    sz = size(RC);
        t = 1:sz(1);
figure(8);

set(gcf,'name','Original time series X and reconstruction RC')
clf;
subplot(2,1,1)
%plot(t,X,'b-',t,sum(RC(:,:),2),'r-');
%legend('Original','Complete reconstruction');
plot(t,real(Xo),'b-',t,real(sum(RC(:,NS:NP),2)),'r-');
legend('Original','Reconstruction with RCs');

subplot(2,1,2)
plot(t,imag(Xo),'b-',t,imag(sum(RC(:,NS:NP),2)),'r-');
legend('Original','Reconstruction with RCs');
saveas(gcf, fullfile(saveDir, sprintf('reconSSA_timeseriesX_reconstructionRCs_%d-%d_%s.png', NS, NP, position)));
end

end