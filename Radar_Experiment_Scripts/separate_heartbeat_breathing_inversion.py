import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import correlate, butter, filtfilt, periodogram

# Load the CSV file
# file_path = 'Data/080724/MC_frbr2_tweaking_RCs/POST_SSA_Data_RCs2-20_Chest.csv'  # Replace with your file path
file_path = 'Data/080724/MC_frbr7/POST_SSA_Data_Chest.csv'  # Replace with your file path
data = pd.read_csv(file_path)

# Extract relevant columns
ecg = data['ECG']
left_radar = data['Left_Radar']
right_radar = data['Right_Radar']
sample_num = data['sample number']

# Effective sampling frequency after decimation
fs = 100  # Adjusted for 10 kHz original sampling rate and 100x decimation

# Scale ECG to have a maximum amplitude of 0.3 and offset it below the radar signals
ecg_scaled = ecg * (0.3 / ecg.abs().max())  # Scale ECG to amplitude of 0.3
min_radar_value = min(left_radar.min(), right_radar.min())
ecg_offset = min_radar_value - 0.15  # Offset ECG below radar signals
ecg_scaled = ecg_scaled + ecg_offset

# Step 1: Cross-Correlation Analysis (Original and Flipped Signals)
def compute_cross_correlation(ecg, radar_signal):
    correlation = correlate(ecg, radar_signal, mode='full')
    lags = np.arange(-len(ecg) + 1, len(ecg))
    return lags, correlation

# Cross-correlation with original and flipped Left and Right Radar signals
lags_left, corr_left = compute_cross_correlation(ecg, left_radar)
lags_left_flipped, corr_left_flipped = compute_cross_correlation(ecg, -left_radar)
lags_right, corr_right = compute_cross_correlation(ecg, right_radar)
lags_right_flipped, corr_right_flipped = compute_cross_correlation(ecg, -right_radar)

# Step 2: Determine whether to use original or flipped signals based on cross-correlation
# Check if the flipped version has a more prominent positive peak
def should_flip_signal(corr, corr_flipped):
    return np.max(corr_flipped) > np.max(corr)

flip_left_radar = should_flip_signal(corr_left, corr_left_flipped)
flip_right_radar = should_flip_signal(corr_right, corr_right_flipped)

# Apply flipping if needed
if flip_left_radar:
    left_radar = -left_radar
if flip_right_radar:
    right_radar = -right_radar

# Step 3: Filtering (Bandpass Filter for Heart Rate Range)
def bandpass_filter(data, lowcut, highcut, fs, order=4):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    y = filtfilt(b, a, data)
    return y

# Define heart rate frequency range (e.g., 1-2 Hz for heartbeats)
lowcut = 0.9
highcut = 2

# Apply bandpass filter to radar signals (with any necessary flipping applied)
left_radar_filtered = bandpass_filter(left_radar, lowcut, highcut, fs)
right_radar_filtered = bandpass_filter(right_radar, lowcut, highcut, fs)

# Step 4: Fourier Transform Analysis
def compute_frequency_spectrum(signal, fs):
    freqs, psd = periodogram(signal, fs)
    return freqs, psd

# Compute frequency spectra for ECG, Left Radar, and Right Radar
freqs_ecg, psd_ecg = compute_frequency_spectrum(ecg, fs)
freqs_left, psd_left = compute_frequency_spectrum(left_radar, fs)
freqs_right, psd_right = compute_frequency_spectrum(right_radar, fs)

# Plot all comparisons in a 2x3 grid
fig, axs = plt.subplots(2, 3, figsize=(18, 10))

# Time-Domain Plot of Signals (with inversion applied if needed)
axs[0, 0].plot(sample_num, left_radar, label='Left Radar', color='blue')
axs[0, 0].plot(sample_num, right_radar, label='Right Radar', color='green')
axs[0, 0].plot(sample_num, ecg_scaled, label='ECG (scaled and offset)', color='red')
axs[0, 0].set_xlabel('Sample Number')
axs[0, 0].set_ylabel('Signal Value')
axs[0, 0].set_title('Time-Domain Signals (Radar and ECG)')
axs[0, 0].legend()

# Cross-Correlation with Left Radar (Original and Flipped)
axs[0, 1].plot(lags_left, corr_left, label='Left Radar', color='blue')
axs[0, 1].plot(lags_left, corr_left_flipped, label='Left Radar (Flipped)', color='orange')
axs[0, 1].set_xlabel('Lag')
axs[0, 1].set_ylabel('Correlation')
axs[0, 1].set_title('Cross-Correlation with Left Radar')
axs[0, 1].legend()

# Cross-Correlation with Right Radar (Original and Flipped)
axs[1, 0].plot(lags_right, corr_right, label='Right Radar', color='green')
axs[1, 0].plot(lags_right, corr_right_flipped, label='Right Radar (Flipped)', color='red')
axs[1, 0].set_xlabel('Lag')
axs[1, 0].set_ylabel('Correlation')
axs[1, 0].set_title('Cross-Correlation with Right Radar')
axs[1, 0].legend()

# Filtered Signals (with inversion applied if needed)
axs[1, 1].plot(sample_num, left_radar_filtered, label='Left Radar (Filtered)', color='blue')
axs[1, 1].plot(sample_num, right_radar_filtered, label='Right Radar (Filtered)', color='green')
axs[1, 1].plot(sample_num, ecg_scaled, label='ECG (scaled and offset)', color='red')
axs[1, 1].set_xlabel('Sample Number')
axs[1, 1].set_ylabel('Signal Value')
axs[1, 1].set_title('Filtered Radar Signals and ECG')
axs[1, 1].legend()

# Frequency Spectrum
axs[0, 2].plot(freqs_ecg, psd_ecg, label='ECG', color='red')
axs[0, 2].plot(freqs_left, psd_left, label='Left Radar', color='blue')
axs[0, 2].plot(freqs_right, psd_right, label='Right Radar', color='green')
axs[0, 2].set_xlabel('Frequency (Hz)')
axs[0, 2].set_ylabel('Power Spectral Density')
axs[0, 2].set_title('Frequency Spectrum of ECG and Radar Signals')
axs[0, 2].set_xlim([0, 5])  # Limit frequency range to focus on heart-related frequencies
axs[0, 2].legend()

# Hide last empty subplot (bottom-right)
axs[1, 2].axis('off')

# Adjust layout and show plot
plt.tight_layout()
plt.show()
