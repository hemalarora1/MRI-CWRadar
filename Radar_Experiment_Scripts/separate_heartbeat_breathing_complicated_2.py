import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import correlate, butter, filtfilt, periodogram, hilbert
from scipy.fft import fft, fftfreq

# Load the CSV file
# file_path = 'Data/080724/MC_frbr2_tweaking_RCs/POST_SSA_Data_RCs2-20_Chest.csv'  # Replace with your file path
# file_path = 'Data/080724/MC_brh1/POST_SSA_Data_Chest.csv'  # Replace with your file path
# file_path = 'Data/081824/KC_frbr4/POST_SSA_Data_Chest.csv'
# file_path = 'Data/081824/KC_frbr1/POST_SSA_Data_Neck.csv'


# Successful filepaths for separating breathing and heartbeat:
file_path = 'Data/081824/KC_frbr1/POST_SSA_Data_Chest.csv'
data = pd.read_csv(file_path)

# Extract relevant columns
ecg = data['ECG']
left_radar = data['Left_Radar']
right_radar = data['Right_Radar']
sample_num = data['sample number']  # Changed to "sample number"

# Effective sampling frequency after decimation
fs = 100  # Adjusted for 10 kHz original sampling rate and 100x decimation

# Scale ECG to have a maximum amplitude of 0.3 and apply a smaller offset below the radar signals
ecg_scaled = ecg * (0.3 / ecg.abs().max())  # Scale ECG to amplitude of 0.3
min_radar_value = min(left_radar.min(), right_radar.min())
ecg_offset = min_radar_value - 0.05  # Offset ECG just below radar signals with a smaller buffer
ecg_scaled = ecg_scaled + ecg_offset

# Function to compute frequency spectrum
def compute_frequency_spectrum(signal, fs):
    freqs, psd = periodogram(signal, fs)
    return freqs, psd

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

# Step 3: Filtering for Cardiac Signal Isolation in Heart and Respiration Bands
def bandpass_filter(data, lowcut, highcut, fs, order=4):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    y = filtfilt(b, a, data)
    return y

# Define heart rate, respiration, and control frequency ranges
heart_rate_band = (1.0, 2.0)
respiration_band = (0.2, 0.5)
control_band = (2.1, 3.0)

# Apply bandpass filters for heart rate, respiration, and control bands
left_radar_heart = bandpass_filter(left_radar, *heart_rate_band, fs)
right_radar_heart = bandpass_filter(right_radar, *heart_rate_band, fs)
# No filtering for ECG since it's the ground truth
ecg_heart = ecg_scaled  

left_radar_resp = bandpass_filter(left_radar, *respiration_band, fs)
right_radar_resp = bandpass_filter(right_radar, *respiration_band, fs)

left_radar_control = bandpass_filter(left_radar, *control_band, fs)
right_radar_control = bandpass_filter(right_radar, *control_band, fs)

# Step 4: Phase Analysis in the Heart Rate Band
def compute_phase_difference(signal1, signal2):
    # Get the analytical signal (complex representation)
    analytic_signal1 = hilbert(signal1)
    analytic_signal2 = hilbert(signal2)
    
    # Calculate phase angles
    phase1 = np.angle(analytic_signal1)
    phase2 = np.angle(analytic_signal2)
    
    # Compute phase difference
    phase_diff = phase1 - phase2
    phase_diff = (phase_diff + np.pi) % (2 * np.pi) - np.pi  # Normalize to [-pi, pi]
    
    return phase_diff

# Compute phase differences between ECG and radar signals in the heart rate band
phase_diff_left = compute_phase_difference(ecg_heart, left_radar_heart)
phase_diff_right = compute_phase_difference(ecg_heart, right_radar_heart)

# Step 5: Signal Quality Assessment (Signal-to-Noise Ratio in Heart Rate Band)
def compute_snr(signal, fs, lowcut, highcut):
    # Perform FFT
    fft_vals = fft(signal)
    fft_freqs = fftfreq(len(signal), 1 / fs)
    
    # Power in the heart rate band
    band_power = np.sum(np.abs(fft_vals[(fft_freqs >= lowcut) & (fft_freqs <= highcut)])**2)
    
    # Total power outside the heart rate band
    noise_power = np.sum(np.abs(fft_vals[(fft_freqs < lowcut) | (fft_freqs > highcut)])**2)
    
    # Calculate SNR
    snr = 10 * np.log10(band_power / noise_power)
    return snr

snr_left = compute_snr(left_radar_heart, fs, *heart_rate_band)
snr_right = compute_snr(right_radar_heart, fs, *heart_rate_band)

# Print SNR values
print(f"SNR (Left Radar, Heart Band): {snr_left:.2f} dB")
print(f"SNR (Right Radar, Heart Band): {snr_right:.2f} dB")

# Plot all comparisons in a 2x3 grid
fig, axs = plt.subplots(2, 3, figsize=(18, 10))

# Time-Domain Plot of Signals
axs[0, 0].plot(sample_num, left_radar, label='Left Radar', color='blue')
axs[0, 0].plot(sample_num, right_radar, label='Right Radar', color='green')
axs[0, 0].plot(sample_num, ecg_scaled, label='ECG (scaled and offset)', color='red')
axs[0, 0].set_xlabel('Sample Number')
axs[0, 0].set_ylabel('Signal Value')
axs[0, 0].set_title('Time-Domain Signals (Radar and ECG)')
axs[0, 0].legend()

# Cross-Correlation with Left Radar
axs[0, 1].plot(lags_left, corr_left, label='Left Radar', color='blue')
axs[0, 1].plot(lags_left, corr_left_flipped, label='Left Radar (Flipped)', color='orange')
axs[0, 1].set_xlabel('Lag')
axs[0, 1].set_ylabel('Correlation')
axs[0, 1].set_title('Cross-Correlation with Left Radar')
axs[0, 1].legend()

# Cross-Correlation with Right Radar
axs[0, 2].plot(lags_right, corr_right, label='Right Radar', color='green')
axs[0, 2].plot(lags_right, corr_right_flipped, label='Right Radar (Flipped)', color='red')
axs[0, 2].set_xlabel('Lag')
axs[0, 2].set_ylabel('Correlation')
axs[0, 2].set_title('Cross-Correlation with Right Radar')
axs[0, 2].legend()

# Filtered Signals in Heart Rate Band (Use scaled ECG as ground truth without additional filtering)
axs[1, 0].plot(sample_num, left_radar_heart, label='Left Radar (Heart Band)', color='blue')
axs[1, 0].plot(sample_num, right_radar_heart, label='Right Radar (Heart Band)', color='green')
axs[1, 0].plot(sample_num, ecg_scaled, label='ECG (scaled and offset)', color='red')  # Scaled ECG as ground truth
axs[1, 0].set_xlabel('Sample Number')
axs[1, 0].set_ylabel('Signal Value')
axs[1, 0].set_title('Filtered Signals (Heart Rate Band)')
axs[1, 0].legend()

# Phase Difference Analysis
axs[1, 1].plot(sample_num, phase_diff_left, label='Phase Diff with Left Radar', color='blue')
axs[1, 1].plot(sample_num, phase_diff_right, label='Phase Diff with Right Radar', color='green')
axs[1, 1].set_xlabel('Sample Number')
axs[1, 1].set_ylabel('Phase Difference (radians)')
axs[1, 1].set_title('Phase Difference in Heart Rate Band')
axs[1, 1].legend()

# Frequency Spectrum in Control Band
freqs_left_control, psd_left_control = compute_frequency_spectrum(left_radar_control, fs)
freqs_right_control, psd_right_control = compute_frequency_spectrum(right_radar_control, fs)
axs[1, 2].plot(freqs_left_control, psd_left_control, label='Left Radar (Control Band)', color='blue')
axs[1, 2].plot(freqs_right_control, psd_right_control, label='Right Radar (Control Band)', color='green')
axs[1, 2].set_xlabel('Frequency (Hz)')
axs[1, 2].set_ylabel('Power Spectral Density')
axs[1, 2].set_title('Frequency Spectrum in Control Band')
axs[1, 2].set_xlim([0, 5])
axs[1, 2].legend()

# Adjust layout to avoid overlaps
plt.tight_layout()
plt.show()