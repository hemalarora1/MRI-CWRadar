import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file
file_path = 'Data/080724/MC_frbr9_vary_RCs/POST_SSA_Data_RCs5-20_Chest.csv'
data = pd.read_csv(file_path)

# Calculate dynamic offset
min_radar_value = min(data['Left_Radar'].min(), data['Right_Radar'].min())
ecg_offset = min_radar_value - (data['ECG'].max() - data['ECG'].min()) * 0.1  # Offset below radar signals

# Scale and shift ECG line to appear below the radar lines dynamically
ecg_scaled = 0.1 * data['ECG'] * 0.5 + ecg_offset  # Scale and apply dynamic offset

# Plotting the data
plt.figure(figsize=(10, 6))

# Plot Left Radar, Right Radar, and scaled ECG
plt.plot(data['sample number'], -1*data['Left_Radar'], label='Left Radar', color='blue')
plt.plot(data['sample number'],  data['Right_Radar'], label='Right Radar', color='green')
plt.plot(data['sample number'], ecg_scaled, label='ECG (scaled and offset)', color='red')

# Adding labels and title
plt.xlabel('Sample Number')
plt.ylabel('Signal Value')
plt.title(file_path)
plt.legend()

# Display the plot
plt.show()