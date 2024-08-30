function Zssa = recMSSA(RC, NS, NP, ShowFigs )
% reconstruct time series from mSSA components NS through NP.
% example:
%    Drecon = recMSSA(RC, 2, 10, 0) ;
%    uses 2nd through 10th mSSA componnents. No figure shown ,
%
%   Drecon = recMSSA(RC, 1, 20, 1) ;
%    uses first 20 mSSA components to reconstruct data. Figures are shown.
%
dimRC = size(RC);
N = dimRC(1) ;   % length of each time series - 
if (length(dimRC)==2)
    nCH =1;
else
    nCH = dimRC(3);  % number of channels
end
M = dimRC(2)/nCH;  % window size for SSA temporal direction

Zssa = zeros(N,nCH);
Xo = Zssa;
for k=1:nCH
    Zssa(:,k) = sum(RC(:,NS:NP,k),2);
    Xo(:,k) = sum(RC(:,:,k),2);
end

if(ShowFigs)
    t = 1:N ;
    figure(8);
    set(gcf,'name','Original time series X and reconstruction RC')
    clf;
    for m=1:nCH
        subplot(2,nCH,m)
        plot(t,real(Xo(:,m)),'b',t,real(Zssa(:,m)),'r'); 
    end
    legend('Orig-real','RCs NS:NP');
    for m=1:nCH
        subplot(2,nCH,m+nCH)
        plot(t,imag(Xo(:,m)),'b',t,imag(Zssa(:,m)),'r');   
    end
    legend('Orig-imag','RCs NS:NP');
end  %if

end  %function