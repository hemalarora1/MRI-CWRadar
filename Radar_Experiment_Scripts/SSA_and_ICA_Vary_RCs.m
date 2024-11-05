% This script handles running post-processing for raw radar data in
% left/right radar pairs, including a full SSA/ICA pipeline with figures
% and processed data saved at the given data filepath.
% Script is adapted from plotGS.m (for SSA) and plotFigs.m (for ICA)

% Script expects data_filepath, I1, I2, Q1, Q2, Radar1, Radar2
% variables to be set in the workspace (from python file).
% If manual inputs are needed:

% data_filepath = "Data/080724/MC_frbr9_vary_RCs/MC_frbr9.csv";
data_filepath = "Data/080724/MC_frbr10_vary_RCs/MC_frbr10.csv";
 
choice = 2;

if (choice == 1)
    % Neck
    I1 = 4;
    Q1 = 5;
    Radar1 = "Left Neck";
    I2 = 6;
    Q2 = 7;
    Radar2 = "Right Neck";
    ECG = 2;
    PPG = 3;
    position = "Neck";
    RC1 = 2;
    RC2 = 20;
else
    % Chest
    I1 = 8;
    Q1 = 9;
    Radar1 = "Left Chest";
    I2 = 10;
    Q2 = 11;
    Radar2 = "Right Chest";
    ECG = 2;
    PPG = 3;
    position = "Chest";
    RC1 = 2;
    RC2 = 20;
end

data = readmatrix(data_filepath, 'NumHeaderLines', 1);  % Ignore the header line
tt = data(:, 1) ;
Ia = data(:, I1); 
Qa = data(:, Q1); 
Ib = data(:, I2); 
Qb = data(:, Q2); 
NEX = 100;
[saveDir, ~, ~] = fileparts(data_filepath); % get directory path for saving figures

L = length(Ia);
PH = 180/pi*angle(Ia+j*Qa);
Za = Ia+j*Qa;
Zb = Ib+j*Qb;

t = zeros(round(L/NEX),1);
n = 1;

% average the data ie boxcar
for ll=1:NEX:L
    t(n) = tt(ll);
    n = n+1;
end

% Decimate Data
XXa = decimate(Za,NEX/10) ;
XXa = decimate(XXa,NEX/10) ;% Chebyshev IIR order 8
XXb = decimate(Zb, NEX/10) ;
XXb = decimate(XXb, NEX/10) ;

PHxa = 180/pi*angle(XXa);
PHxb = 180/pi*angle(XXb);
Ixa = real(XXa); Qxa = imag(XXa);
Ixb = real(XXb); Qxb = imag(XXb);

% X = X-mean(X);
N = length(XXa);
M = round(N/4);  % just window size

% perform SSA on data to clean it up
[RCa, LAMBDA] = compSSA(XXa, 4, 1, saveDir, Radar1);
[RCb, LAMBDA] = compSSA(XXb, 4, 1, saveDir, Radar2);
% reconstructs data from components 2 through 20
Arec = reconSSA(RCa, RC1, RC2, 1, saveDir, Radar1);
Brec = reconSSA(RCb, RC1, RC2, 1, saveDir, Radar2);

% do PCA + rotation on cleaned up SSA data
dat_a = Arec;
dat_b = Brec;
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
Vca(:,1) = -Vca(:,1) ; % convert first vector to opposite direction?
Vca = -Vca;
Vcb(:,1) = -Vcb(:,1) ;
dprin_a = (Vca'*[real(dtm_a), imag(dtm_a)]')';
dprin_b = (Vcb'*[real(dtm_b), imag(dtm_b)]')';

% Now attempt to fit a circle of radius Ro onto the data.
% here brute force approxmate shift to find origin, that minimizes
% variance of fit. Pretty clumsy approach.
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

% ECG Data
% WARNING - I is PPG, Q is ECG

I = data(:, PPG);
Q = data(:, ECG);

L = length(I);
PH = 180/pi*angle(I+j*Q);
Z = I+j*Q;

t = zeros(round(L/NEX),1);
n = 1;

%average the data ie boxcar
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

figure(9);
set(gcf,'name','Time series versus ECG Data');
clf;

% PLOT SSA Result Against ECG baseline
% Push ecg line below other signals
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
saveas(gcf, fullfile(saveDir, sprintf('final_plotGS_timeseriesvsECG_RCs%d-%d_%s.png', RC1, RC2, position)));
end

% SAVE SSA DATA
if(1)
% Define the column headers
headers = {'sample number', 'ECG', 'PPG', 'Left_Radar', 'Right_Radar'};

% Specify the filename and path
csv_name = sprintf("%s/POST_SSA_Data_RCs%d-%d_%s.csv", saveDir, RC1, RC2, position);

% Write the headers and then append the data
fid = fopen(csv_name, 'w');  % Open file for writing
fprintf(fid, '%s,', headers{1:end-1});  % Write headers except the last one
fprintf(fid, '%s\n', headers{end});  % Write the last header followed by a newline
fclose(fid);

% Now append the data
samples = (0:size(Qx, 1) - 1)';
post_SSA_data = [samples, Qx, Ix, dprin_a(:, 2), dprin_b(:,2)];
writematrix(post_SSA_data, csv_name, 'WriteMode', 'append');
end
% -------------------------------------------------------------------------
% FASTICA/RICA START
% -------------------------------------------------------------------------

% script needs access to fastica function - keep fastica_25 folder in
% same directory as this script
addpath('./FastICA_25');

% post-SSA data
chA = dprin_a(:,2);
chB = dprin_b(:,2);

comb_radars = [chA,chB];
fica_radars = fastica(comb_radars'); fica_radars = fica_radars'; 

wdat = prewhiten(comb_radars);

% important to know if source is super (+1) or sub(-1) gaussian
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

rica_model = rica(wdat, 2, 'NonGaussianityIndicator', nonGaussianityIndicator);
rica_radars = transform(rica_model, wdat);
fica_radars_prewhiten = fastica(wdat'); fica_radars_prewhiten = fica_radars_prewhiten';

if(1)
figure(1);
set(gcf, 'Position', [100, 100, 1400, 1000]);  % [left, bottom, width, height]

subplot(5,1,1);
plot(comb_radars);
title("Post-SSA");
legend(Radar1, Radar2);

subplot(5,1,2);
plot(rica_radars);
title("RICA");
legend("IC1", "IC2");

subplot(5,1,3);
plot(fica_radars);
title("FASTICA");
legend("IC1", "IC2");

subplot(5,1,4);
plot(fica_radars_prewhiten);
title("FASTICA-PREWHITENED");
legend("IC1", "IC2")

subplot(5,1,5);
plot(Qx);
title("ECG");
legend("ECG")
saveas(gcf, fullfile(saveDir, sprintf('RICA_FASTICA_vs_ECG_RCs%d-%d_%s.png', RC1, RC2, position)));
end

% Define the filename with the position variable
csv_name_ica = sprintf("%s/POST_ICA_DATA_RCs%d-%d_%s.csv", saveDir, RC1, RC2, position);

% Define column headers for the FastICA data
headers_ica = {'sample_number', 'IC1_RICA', 'IC2_RICA', 'IC1_FASTICA', 'IC2_FASTICA', 'IC1_FASTICA_PREWHITENED', 'IC2_FASTICA_PREWHITENED'};

% Write the headers to the CSV file
fid_ica = fopen(csv_name_ica, 'w');
fprintf(fid_ica, '%s,', headers_ica{1:end-1});
fprintf(fid_ica, '%s\n', headers_ica{end});
fclose(fid_ica);

% Create a sample index
samples_ica = (0:size(rica_radars, 1) - 1)';

% Combine FastICA and RICA data into a matrix
post_ica_data = [samples_ica, rica_radars, fica_radars, fica_radars_prewhiten];

% Append the data to the CSV file
writematrix(post_ica_data, csv_name_ica, 'WriteMode', 'append');