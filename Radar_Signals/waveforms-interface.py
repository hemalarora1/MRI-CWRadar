"""
    scope 1 (ADP3450): neck radars
        left neck on channel 1 and 2
        right neck on channel 3 and 4
    scope2 (ADP3450): chest radars
        left chest on channel 1 and 2
        right chest on channel 3 and 4
    scope3 (AD2):
        ecg on channel 1
        ppg on channel 2  
        
    trigger 1 of scope 1 connected to trigger 1 of scope 2
    trigger 2 of scope 2 connected to trigger 1 of scope 3
    
    waveforms must be installed on laptop to run script

    Author: Hemal Arora (REU Summer 2024) - hemal1@stanford.edu
"""

from ctypes import *
from dwfconstants import *
import time
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

# device serial numbers
SCOPE1_SN = "210018B04928"
SCOPE2_SN = "210018B9D8A8"
# SCOPE3_SN = "" #AD2
SUPPLY_SN = "210395B7C440"

# experiment parameters
SAMPLING_FREQUENCY = 100000.0 # 100 kHz sampling rate
RECORD_LENGTH = 10 # 10 second recording length
NUM_SAMPLES = RECORD_LENGTH * SAMPLING_FREQUENCY
OFFSET = -2.0 # -2V offset
CHANNEL_RANGE = 5.0 # +/- 5V input range
V_PLUS = 5.0 # +5V Supply
V_MINUS = -5.0 # -5V Supply
CURRENT_LIMIT = 0.4 # 10 mA to 3A range available - check with Dr. Scott

# load waveforms library
current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
dwf = cdll.LoadLibrary(f"{current_directory}/Waveforms-SDK/dwf")

class device:
    def __init__(self, id, serial_number, description):
        self.serial_number = serial_number
        self.id = id
        self.handle = c_int()
        self.description = description
        self.acqstatus = c_byte(0)

def open_devices():
    # Initialize device variables
    scope1 = None
    scope2 = None
    supply = None
    device_num = c_int(0)

    # enumerate all connected devices
    dwf.FDwfEnum(c_int(0), byref(device_num))
    print("Number of Devices: " + str(device_num.value))

    # generate device objects with assigned enumeration IDs
    current_sn = create_string_buffer(16)

    for device_id in range(device_num.value):
        dwf.FDwfEnumSN(c_int(device_id), current_sn)
        sn = current_sn.value.decode('utf-8') # convert to python string for comparison
        sn = sn[3:] # remove "SN:" prefix
        
        if sn == SCOPE1_SN:
            scope1 = device(id=device_id, serial_number=SCOPE1_SN, description="ADP3450 Neck Radar") # neck radar scope
        elif sn == SCOPE2_SN:
            scope2 = device(id=device_id, serial_number=SCOPE2_SN, description="ADP3450 Chest Radar") # chest radar scope
        # if sn == SCOPE3_SN:
        #     scope1 = device(id=device_id, serial_number=SCOPE3_SN, description="AD2 ECG/PPG Sensor") # ECG/PPG scope
        elif sn == SUPPLY_SN:
            supply = device(id=device_id, serial_number=SUPPLY_SN, description="DPS3340 Power Supply") # power supply
        else:
            raise Exception(f"Unknown Device Serial Number Found: {sn}")

    # opening all devices
    for component in [scope1, scope2, supply]: # replace with [scope1, scope2, scope3, supply] after setting up AD2
        dwf.FDwfDeviceOpen(c_int(component.id), byref(component.handle))
        if component.handle.value == hdwfNone.value:
            szerr = create_string_buffer(512)
            dwf.FDwfGetLastErrorMsg(szerr)
            print(szerr.value)
            print("Failed to open device!")
            quit()
        dwf.FDwfDeviceAutoConfigureSet(component.handle, c_int(0)) # 0 = the device will only be configured when FDwf###Configure is called
        print(f"{component.description} - Opened Successfully!")

    return scope1, scope2, supply

def power_supply_config(supply):
    # configuring power supply
    dwf.FDwfAnalogIOChannelNodeSet(supply.handle, c_int(1), c_int(0), c_double(True)) # enable output 2
    dwf.FDwfAnalogIOChannelNodeSet(supply.handle, c_int(1), c_int(1), c_double(V_MINUS)) # set v_minus
    dwf.FDwfAnalogIOChannelNodeSet(supply.handle, c_int(1), c_int(2), c_double(-1 * CURRENT_LIMIT)) 
    dwf.FDwfAnalogIOChannelNodeSet(supply.handle, c_int(2), c_int(0), c_double(True)) # enable output 3
    dwf.FDwfAnalogIOChannelNodeSet(supply.handle, c_int(2), c_int(1), c_double(V_PLUS)) # set v_plus
    dwf.FDwfAnalogIOChannelNodeSet(supply.handle, c_int(2), c_int(2), c_double(CURRENT_LIMIT))
    dwf.FDwfAnalogIOEnableSet(supply.handle, c_int(True)) # master enable on
    dwf.FDwfAnalogIOConfigure(supply.handle)

def adp_config(scope1, scope2):
    # configuring ADP3450s
    for adp in [scope1, scope2]:
        # acquisition settings for analog in channels
        dwf.FDwfAnalogInAcquisitionModeSet(adp.handle, acqmodeRecord)
        dwf.FDwfAnalogInChannelEnableSet(adp.handle, c_int(-1), c_int(1))
        dwf.FDwfAnalogInChannelRangeSet(adp.handle, c_int(-1), c_double(CHANNEL_RANGE)) # set channel range for all channels
        dwf.FDwfAnalogInChannelOffsetSet(adp.handle, c_int(-1), c_int(OFFSET)) # -2V offset on all channels
        dwf.FDwfAnalogInFrequencySet(adp.handle, c_double(SAMPLING_FREQUENCY))
        dwf.FDwfAnalogInRecordLengthSet(adp.handle, c_double(RECORD_LENGTH))
        dwf.FDwfAnalogInBufferSizeSet(adp.handle, c_int(NUM_SAMPLES))

        # trigger settings - if this doesn't work use built in trigsrcNone
        if adp.serial_number == SCOPE1_SN: # scope 1 used as master device
            # on first device drive External T1 (0) with trigsrcPC
            dwf.FDwfAnalogInTriggerSourceSet(adp.handle, trigsrcPC) # second device is triggered based on pin (may need trigsrcPC for device 1)
            dwf.FDwfDeviceTriggerSet(adp.handle, c_int(0), trigsrcAnalogIn) # I think this sets trigger 1 as an output based on trigsrcPC?
        else:
            dwf.FDwfAnalogInTriggerSourceSet(adp.handle, trigsrcExternal1) # second device is triggered based on pin (may need trigsrcPC for device 1)
            dwf.FDwfDeviceTriggerSet(adp.handle, c_int(1), trigsrcAnalogIn)
        dwf.FDwfAnalogInTriggerAutoTimeoutSet(adp.handle, c_double(0)) # disable auto trigger
        dwf.FDwfAnalogInTriggerPositionSet(adp.handle, c_double(0.0))
        dwf.FDwfAnalogInConfigure(adp.handle, c_int(True), c_int(False)) # configure (maybe don't reset auto-trigger?) and don't start reading data?
        print(f"{adp.description}: Configured!")
        time.sleep(3) # wait for offset to stabilise
        
def ad2_config(scope3):
    # acquisition settings for analog in channels
    dwf.FDwfAnalogInAcquisitionModeSet(scope3.handle, acqmodeRecord)
    dwf.FDwfAnalogInChannelEnableSet(scope3.handle, c_int(-1), c_int(1))
    dwf.FDwfAnalogInChannelRangeSet(scope3.handle, c_int(-1), c_double(CHANNEL_RANGE)) # set channel range for all channels
    dwf.FDwfAnalogInFrequencySet(scope3.handle, c_double(SAMPLING_FREQUENCY))
    dwf.FDwfAnalogInRecordLengthSet(scope3.handle, c_double(RECORD_LENGTH))
    dwf.FDwfAnalogInBufferSizeSet(scope3.handle, c_int(NUM_SAMPLES))
    
    # trigger settings
    dwf.FDwfAnalogInTriggerSourceSet(scope3.handle, trigsrcExternal1)
    dwf.FDwfAnalogInTriggerAutoTimeoutSet(scope3.handle, c_double(0)) #disable auto trigger
    dwf.FDwfAnalogInTriggerPositionSet(scope3.handle, c_double(0.0))
    dwf.FDwfAnalogInConfigure(scope3.handle, c_int(True), c_int(False)) # configure (maybe don't reset auto-trigger?) and don't start reading data
    print(f"{scope3.description}: Configured!")
    time.sleep(3)

def extract_scope_data(handle):
    data = []
    channel_num = c_int(0)
    samples_buffer = (c_double() * NUM_SAMPLES)()
    dwf.FDwfAnalogInChannelCount(handle, byref(channel_num))
    for channel in range(channel_num.value):
        dwf.FDwfAnalogInStatusData(handle, channel, samples_buffer, NUM_SAMPLES)
        data.append(list(samples_buffer))
    return data

def run_experiment(scope1, scope2):
    dwf.FDwfDeviceTriggerPC(scope1.handle) # trigger master device (which triggers secondary devices)
    
    while True: # recording in progress - do I need to track corrupted/lost data?
        dwf.FDwfAnalogInStatus(scope1.handle, 1, byref(scope1.acqstatus))
        dwf.FDwfAnalogInStatus(scope2.handle, 1, byref(scope2.acqstatus))
        # dwf.FDwfAnalogInStatus(scope3.handle, 1, byref(status3)) # for ad2
        if scope1.acqstatus.value == DwfStateDone.value and scope2.acqstatus.value == DwfStateDone.value:
            print("Acquisition done!")
            break
        time.sleep(0.5)
    
    neck_radar_data = extract_scope_data(scope1.handle)
    chest_radar_data = extract_scope_data(scope2.handle)
    # ecg_data = extract_scope_data(scope3.handle)
    
    experiment_data = {}
    # experiment_data["ECG"] = ecg_data[0]
    # experiment_data["PPG"] = ecg_data[1]
    experiment_data["Left Neck - I"] = neck_radar_data[0]
    experiment_data["Left Neck - Q"] = neck_radar_data[1]
    experiment_data["Right Neck - I"] = neck_radar_data[2]
    experiment_data["Right Neck - Q"] = neck_radar_data[3]
    experiment_data["Left Chest - I"] = chest_radar_data[0]
    experiment_data["Left Chest - Q"] = chest_radar_data[1]
    experiment_data["Right Chest - I"] = chest_radar_data[2]
    experiment_data["Right Chest - Q"] = chest_radar_data[3]
    
    final_data = pd.DataFrame(experiment_data)
    
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))
    
    # Plot Neck data
    axs[0].plot(time, experiment_data["Left Neck - I"], label="Left Neck - I")
    axs[0].plot(time, experiment_data["Left Neck - Q"], label="Left Neck - Q")
    axs[0].plot(time, experiment_data["Right Neck - I"], label="Right Neck - I")
    axs[0].plot(time, experiment_data["Right Neck - Q"], label="Right Neck - Q")
    axs[0].set_title("Neck Signals")
    axs[0].legend()

    # Plot Chest data
    axs[1].plot(time, experiment_data["Left Chest - I"], label="Left Chest - I")
    axs[1].plot(time, experiment_data["Left Chest - Q"], label="Left Chest - Q")
    axs[1].plot(time, experiment_data["Right Chest - I"], label="Right Chest - I")
    axs[1].plot(time, experiment_data["Right Chest - Q"], label="Right Chest - Q")
    axs[1].set_title("Chest Signals")
    axs[1].legend()

    # # Plot ECG/PPG data
    # axs[2].plot(time, experiment_data["ECG"], label="ECG")
    # axs[2].plot(time, experiment_data["PPG"], label="PPG")
    # axs[2].set_title("ECG/PPG Signals")
    # axs[2].legend()
    
    # Adjust layout and display
    plt.tight_layout()
    plt.show()
    
    return final_data

def diagnostics(scope1, scope2):
    max_buffer_size = c_int()
    dwf.FDwfAnalogInBufferSizeInfo(scope1.handle, None, byref(max_buffer_size))
    print(f"Max buffer size (scope1): {max_buffer_size.value}")

    current_buffer_size = c_int()
    dwf.FDwfAnalogInBufferSizeGet(scope1.handle, byref(current_buffer_size))
    print(f"Current buffer size (scope 1): {current_buffer_size.value}")

    trigger_position = c_double()
    dwf.FDwfAnalogInTriggerPositionGet(scope1.handle, byref(trigger_position))
    print(f"Current trigger position (scope1): {trigger_position.value}")
    
def main():
    scope1, scope2, supply = open_devices()
    power_supply_config(supply)
    adp_config(scope1, scope2)
    diagnostics(scope1, scope2)
    
    while True:
        print("\n\nMenu:")
        print("1. Run Experiment")
        print("2. Quit")
        choice = input("Enter your choice: ")
        
        while choice != '1' and choice != '2':
            print("Invalid Choice, Try Again...")
            choice = input("Enter your choice: ")
        
        if choice == '1':
            experiment_data = run_experiment()
            save = input("Do you want to save the experiment data (y/n)?: ")
            while save != 'y' and choice != 'n':
                print("Invalid Choice, Try Again...")
                save = input("Do you want to save the experiment data (y/n)?: ")
            if save == 'y':
                filename = input("Enter the filename to save the data: ")
                experiment_data.to_csv(f"{current_directory}/Data/{filename}.csv", index=False)
                print("Data saved in 'Data' folder")
            if save == 'n':
                continue
        elif choice == '2':
            dwf.FDwfDeviceCloseAll() # close all devices
            break
    
    print("Session Terminated!")

if __name__ == '__main__':
    main()