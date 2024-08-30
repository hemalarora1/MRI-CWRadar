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
[RCa, LAMBDA] = compSSA(XXa, 4, 1, saveDir, Radar1);
[RCb, LAMBDA] = compSSA(XXb, 4, 1, saveDir, Radar2);

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


%PLOTS

figure(9);
set(gcf,'name','Time series versus ECG Data');
clf;

% Try to push ecg line below other signals
% Calculate the minimum values of the signals
ecg = (0.3/(max(Qx) - min(Qx))) * Qx;
max_ecg = max(ecg);
min_left = min(dprin_a(:, 2));
min_right = min(dprin_b(:, 2));

% Find the minimum of the left and right neck/chest signals
min_other_signals = min(min_left, min_right);

% Calculate the offset needed to make the ECG signal maximum below the other signals' minimum
offset = max_ecg - min_other_signals;  % The constant can be adjusted

% Apply the offset to the ECG signal
adjusted_ecg = ecg - offset;

if(1)
plot(adjusted_ecg, 'g');
hold on;
%plot(-dprin_a(:, 2),'b'); % for chest - set a
plot(dprin_a(:, 2),'b');  % for neck -setb
plot(dprin_b(:, 2),'r');
hold off;
legend('ECG', Radar1, Radar2);
saveas(gcf, fullfile(saveDir, sprintf('final_plotGS_timeseriesvsECG_%s.png', position)));
end

% USED FOR PPG/RIGHT NECK COMPARISON ON 080724/MC_BRH1
% figure(10)
% set(gcf, 'Position', [100, 100, 1200, 800]);  % [left, bottom, width, height]
% plot(Ix, 'g');
% hold on;
% plot(dprin_b(:, 2) + 2.6,'r');
% hold off;
% legend('PPG', Radar2);
% saveas(gcf, fullfile(saveDir, sprintf('PPG_vs_right_neck_%s.png', Radar2)));

% SAVE DATA
% 
% % Define the column headers
% headers = {'sample number', 'PPG', 'Right Neck'};
% 
% % Specify the filename and path
% csv_name = sprintf("%s/POST_SSA_PPG_RIGHT_NECK.csv", saveDir)
% 
% % Write the headers and then append the data
% fid = fopen(filename, 'w');  % Open file for writing
% fprintf(fid, '%s,', headers{1:end-1});  % Write headers except the last one
% fprintf(fid, '%s\n', headers{end});  % Write the last header followed by a newline
% fclose(fid);
% 
% % Now append the data
% samples = (0:size(data, 1) - 1)';
% ppg_rn_data = [samples, Ix, dprin_b(:, 2)]
% writematrix(data, filename, 'WriteMode', 'append');
