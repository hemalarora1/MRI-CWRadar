% Data "sampled" at 10KHz
%fast ICA trials for multiple channels

% if manual inputs needed:
% data_filepath = "Data/080724/MC_brh1/MC_brh1.csv";
% 
% choice = 1;
% 
% if (choice == 1)
%     % Neck
%     I1 = 4;
%     Q1 = 5;
%     Radar1 = "Left Neck";
%     I2 = 6;
%     Q2 = 7;
%     Radar2 = "Right Neck";
%     ECG = 2;
%     PPG = 3;
%     position = "Neck";
% else
%     % Chest
%     I1 = 8;
%     Q1 = 9;
%     Radar1 = "Left Chest";
%     I2 = 10;
%     Q2 = 11;
%     Radar2 = "Right Chest";
%     ECG = 2;
%     PPG = 3;
%     position = "Chest";
% end


% Data "sampled" at 10KHz
%PLOTS time series + reconstructed w/SSA w/2 channels
% This script expects data_filepath, I1, I2, Q1, Q2, Radar1, Radar2 variables to be set in the workspace (from python file)

% if manual inputs needed:
% data_filepath = "Data/080724/MC_brh1/MC_brh1.csv";
% 
% choice = 1;
% 
% if (choice == 1)
%     % Neck
%     I1 = 4;
%     Q1 = 5;
%     Radar1 = "Left Neck";
%     I2 = 6;
%     Q2 = 7;
%     Radar2 = "Right Neck";
%     ECG = 2;
%     PPG = 3;
%     position = "Neck";
% else
%     % Chest
%     I1 = 8;
%     Q1 = 9;
%     Radar1 = "Left Chest";
%     I2 = 10;
%     Q2 = 11;
%     Radar2 = "Right Chest";
%     ECG = 2;
%     PPG = 3;
%     position = "Chest";
% end

data = readmatrix(data_filepath, 'NumHeaderLines', 1);  % Ignore the first header line
tt = data(:, 1) ;
Ia = data(:, I1); 
Qa = data(:, Q1); 
Ib = data(:, I2); 
Qb = data(:, Q2); 
NEX = 100;
[saveDir, ~, ~] = fileparts(data_filepath); % get directory path for saving figures

%COMPARED DATA STUFF

%time series stuff (dopproc2)
L = length(Ia);
PH = 180/pi*angle(Ia+j*Qa);
Za = Ia+j*Qa;
Zb = Ib+j*Qb;

t = zeros(round(L/NEX),1);
n = 1;
%average the AD2 data ie boxcar
for ll=1:NEX:L
    t(n) = tt(ll);
    n = n+1;
end

%Decimate Data
XXa = decimate(Za,NEX/10) ;
XXa = decimate(XXa,NEX/10) ;% Chebyshev IIR order 8
XXb = decimate(Zb, NEX/10) ;
XXb = decimate(XXb, NEX/10) ;

PHxa = 180/pi*angle(XXa);
PHxb = 180/pi*angle(XXb);
Ixa = real(XXa); Qxa = imag(XXa);
Ixb = real(XXb); Qxb = imag(XXb);

%X = X-mean(X);
N = length(XXa);
M = round(N/4);  % just window size

%perform SSA on data to clean it up
[RCa, LAMBDA] = compSSA(XXa, 4, 0, saveDir, Radar1);
[RCb, LAMBDA] = compSSA(XXb, 4, 0, saveDir, Radar2);

%FBa
%Arec = reconSSA( RCa, 1, 2, 0) + reconSSA( RCa, 5, 20, 0);
%Brec = reconSSA( RCb, 1, 2, 0) + reconSSA( RCb, 5, 20, 0);

%DBa
%Arec = reconSSA( RCa, 1, 20, 0) + reconSSA( RCa, 3, 20, 0);
%Brec = reconSSA( RCb, 1, 20, 0);

%BH
%Arec = reconSSA( RCa, 2, 20, 0);
%Brec = reconSSA( RCb, 2, 20, 0);

%b
Arec = reconSSA( RCa, 2, 20, 0, saveDir, Radar1) ;     % reconstructs data from components 1 thru 20
Brec = reconSSA( RCb, 2, 20, 0, saveDir, Radar2) ;  % change both to 1,20 only for FB2a

% free breath chest FB2a.txt :use 1:20
% free breath neck FB2b.txt  :use 2:20
% breath hold chest BH2a.txt : use 2:20
% breath hold neck BH2b.txt  :use 2:20


%do PCA + rotation stuff on cleaned up SSA data
dat_a = Arec;
dat_b = Brec; % change to your data set
%dat = dat(2:2098,:);
dtm_a = dat_a-mean(dat_a);  % strips off average
dtm_b = dat_b-mean(dat_b);


% Now find the principle axis for the data
Ca = cov(real(dtm_a), imag(dtm_a));
Cb = cov(real(dtm_b), imag(dtm_b));
%C(1,1) = mean(real(dtm).^2); C(2,2) = mean(imag(dtm).^2); %same result
%C(2,1) = mean(real(dtm).*imag(dtm)); C(1,2) = C(2,1);
[Vca,D] = eig(Ca);  % why is det(Vc) = -1? means mirror image
[Vcb,D] = eig(Cb);
%[U,S,V] = svd(C);
% next line can comment out if mirror image issue doesn't show up.
Vca(:,1)=-Vca(:,1) ; % convert first vector to opposite direction?
Vca = -Vca;
Vcb(:,1)=-Vcb(:,1) ;
dprin_a = ( Vca'*[real(dtm_a), imag(dtm_a)]' )';
dprin_b = ( Vcb'*[real(dtm_b), imag(dtm_b)]' )';

% Now attempt to fit a circle of radius Ro onto the data.
% here brute force approxmate shift to find origin, that minimizes
% variance of fit. Pretty clumsy.
arc_a = dprin_a(:,1)+i*dprin_a(:,2);
arc_b = dprin_b(:,1)+i*dprin_b(:,2);
Ro_a = abs(arc_a);
Ro_b = abs(arc_b);
Fo_a = mean(Ro_a.^2)-mean(Ro_a)^2;
Fo_b = mean(Ro_b.^2)-mean(Ro_b)^2;
ko_a=0;
ko_b=0;
for k_a=-.5:.05:.5
    Ri_a = abs(arc_a+k_a);
    F_a = mean(Ri_a.^2)-mean(Ri_a)^2;
    if (F_a<Fo_a)
        ko_a =k_a;
        Fo_a =F_a;
    end
end
for k_b =-.5:.05:.5
    Ri_b = abs(arc_b+k_b);
    F_b = mean(Ri_b.^2)-mean(Ri_b)^2;
    if (F_b<Fo_b)
        ko_b = k_b;
        Fo_b =F_b;
    end
end



%ECG DATA STUFF
% WARNING - I is PPG, Q is ECG

% load fullFB2c.txt; tt=fullFB2c(:,1); I = fullFB2c(:,2); Q = fullFB2c(:,3); NEX = 100;
%load fullDB2c.txt; tt=fullDB2c(:,1); I = fullDB2c(:,2); Q = fullDB2c(:,3); NEX = 100;
I = data(:, PPG);
Q = data(:, ECG);

L = length(I);
PH = 180/pi*angle(I+j*Q);
Z = I+j*Q;

t = zeros(round(L/NEX),1);
n = 1;
%average the AD2 data ie boxcar
for ll=1:NEX:L
    t(n) = tt(ll);
    n = n+1;
end

%decimate data
XX = decimate(Z,NEX/10) ; % Chebyshev IIR order 8
XX = decimate(XX,NEX/10);
PHx = 180/pi*angle(XX);
Ix = real(XX); Qx = imag(XX);

%X = X-mean(X);
N = length(XX);
M = round(N/4);  % just window size

% -------------------------------------------------------------------------
% DUAL FASTICA START (not properly adapted since it should take output of
% SSA rather than raw data)
% Also outputs are not plotted correctly since they should be labelled as
% independent components
% -------------------------------------------------------------------------

% script needs access to fastica function - keep fastica_25 folder in
% same directory as this script
addpath('./FastICA_25');

L = length(Ia);
PH = 180/pi*angle(Ia+j*Qa);
Za = Ia+j*Qa;
Zb = Ib+j*Qb;

%decimate data
XXa = decimate(Za,NEX/10) ;
XXa = decimate(XXa,NEX/10) ;% Chebyshev IIR order 8
XXb = decimate(Zb, NEX/10) ;
XXb = decimate(XXb, NEX/10) ;

Ixa = real(XXa); Qxa = imag(XXa);
Ixb = real(XXb); Qxb = imag(XXb);

dat_a = XXa;
dat_b = XXb; % change to your data set
%dat = dat(2:2098,:);
dtm_a = dat_a-mean(dat_a);  % strips off average
dtm_b = dat_b-mean(dat_b);


% Now find the principle axis for the data
Ca = cov(real(dtm_a), imag(dtm_a));
Cb = cov(real(dtm_b), imag(dtm_b));
[Vca,D] = eig(Ca);  % why is det(Vc) = -1? means mirror image
[Vcb,D] = eig(Cb);
%[U,S,V] = svd(C);
% next line can comment out if mirror image issue doesn't show up.
% check with dr. scott if I should leave the flipping as is
Vca(:,1)=-Vca(:,1) ; % convert first vector to opposite direction?
Vcb(:,1)=-Vcb(:,1) ;
dprin_a = ( Vca'*[real(dtm_a), imag(dtm_a)]' )';
dprin_b = ( Vcb'*[real(dtm_b), imag(dtm_b)]' )';

%??
t = [0:0.001:1]';
Npts = length(t);

%ordat

%data inputs
%data = [Ixa, Qxa] ; % can we do ICA on a single I/Q sensor?
%data = [dprin_a, dprin_b];  % do ICA on separate principle components?
data = [Ixa, Qxa, Ixb, Qxb]; % ICA on ALL I/Q?
q = 2 ; % number of PCA components
%q = 4;
[coeff,Dpca,latent,tsquared,explained,mu] = pca(data, 'NumComponents',q);

%kurtosis(ordat)    % pure infinite length gaussian is 3
wdat = prewhiten(data);

% hugely important to know if source is super (+1) or sub(-1) gaussian
% dynamically decide if source is super or sub-gaussian by kurtosis
kurtosis_values = kurtosis(wdat);
nonGaussianityIndicator = zeros(1, size(wdat, 2)); % Preallocate
for i = 1:length(kurtosis_values)
    if kurtosis_values(i) > 3
        nonGaussianityIndicator(i) = 1;  % Super-Gaussian
    else
        nonGaussianityIndicator(i) = -1; % Sub-Gaussian
    end
end

%Noisy
Mdl = rica(wdat, 4, 'NonGaussianityIndicator', nonGaussianityIndicator);
%RICA
Dica = transform(Mdl,wdat);
%FASTICA
Fica = fastica(wdat'); Fica = Fica';


% ECG/PPG data
I = raw_data(:, PPG);
Q = raw_data(:, ECG);

L = length(I);
PH = 180/pi*angle(I+j*Q);
Z = I+j*Q;

t = zeros(round(L/NEX),1);
n = 1;
%average the AD2 data ie boxcar
for ll=1:NEX:L
    t(n) = tt(ll);
    n = n+1;
end

%decimate data
XX = decimate(Z,NEX/10) ; % Chebyshev IIR order 8
XX = decimate(XX,NEX/10);
PHx = 180/pi*angle(XX);
Ix = real(XX); Qx = imag(XX);

%X = X-mean(X);
N = length(XX);
M = round(N/4);  % just window size

figure(1);
set(gcf, 'Position', [100, 100, 1200, 800]);  % [left, bottom, width, height]
subplot(2,1,1);
plot(data);
title("Data");
legend(sprintf('%s - I', Radar1), sprintf('%s - Q', Radar1), sprintf('%s - I', Radar2), sprintf('%s - Q', Radar2));
subplot(2,1,2);
plot(wdat);
title("Whitened Data");
legend(sprintf('%s - I', Radar1), sprintf('%s - Q', Radar1), sprintf('%s - I', Radar2), sprintf('%s - Q', Radar2));
saveas(gcf, fullfile(saveDir, sprintf('ICA_original_vs_whitened_%s.png', position)));

figure(2);
set(gcf, 'Position', [100, 100, 1200, 800]);  % [left, bottom, width, height]
%subplot(4,1,1);
%plot(ordat);
%title("Original");
subplot(4,1,1);
plot(data);
title("Noisy");
legend(sprintf('%s - I', Radar1), sprintf('%s - Q', Radar1), sprintf('%s - I', Radar2), sprintf('%s - Q', Radar2));
subplot(4,1,2);
plot(Dica);
title("RICA");
legend(sprintf('%s - I', Radar1), sprintf('%s - Q', Radar1), sprintf('%s - I', Radar2), sprintf('%s - Q', Radar2));
subplot(4,1,3);
plot(Fica);
title("FASTICA");
legend(sprintf('%s - I', Radar1), sprintf('%s - Q', Radar1), sprintf('%s - I', Radar2), sprintf('%s - Q', Radar2));
subplot(4,1,4);
plot(Qx);
title("ECG");
legend("ECG")
saveas(gcf, fullfile(saveDir, sprintf('RICA_FASTICA_vs_ECG_%s.png', position)));

%w/q = 2
%figure(3);
%subplot(2,1,1);
%histogram(data(:,1));
%title("wave a");
%subplot(2,1,2);
%histogram(data(:,2));
%title("wave b");

%w/q = 4
figure(3);
set(gcf, 'Position', [100, 100, 1200, 800]);  % [left, bottom, width, height]
subplot(4,1,1);
histogram(data(:,1));
title(sprintf("wave a: %s - I", Radar1));
subplot(4,1,2);
histogram(data(:,2));
title(sprintf("wave b: %s - Q", Radar1));
subplot(4,1,3);
histogram(data(:,3));
title(sprintf("wave c: %s - I", Radar2));
subplot(4,1,4);
histogram(data(:,4));
title(sprintf("wave d: %s - Q", Radar2));
saveas(gcf, fullfile(saveDir, sprintf('ICA_waves_histograms_%s.png', position)));

