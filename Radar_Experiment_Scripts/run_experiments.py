"""
    Scope 1 (ADP3450): neck radars
        left neck on channel 1 and 2
        right neck on channel 3 and 4
    Scope2 (ADP3450): chest radars
        left chest on channel 1 and 2
        right chest on channel 3 and 4
    Scope3 (AD2):
        ecg on channel 1
        ppg on channel 2  
    
    Trigger 1 of all devices are connected for clock frequency
    Trigger 2 of all devices are connected for trigger signal
    
    Waveforms and Matlab applications must be installed on laptop to run script.
    Waveforms-SDK folder should also remain in this script's directory.
    Waveforms SDK Reference Manual can be used for more info on waveforms SDK functions
    SSA_and_ICA.m is typically used for post-processing and should be in the script directory
    along with compSSA, reconSSA, FastICA_25 folder.

    Author: Hemal Arora (REU Summer 2024) - hemal1@stanford.edu
    Contact me for any questions!
"""

from ctypes import *
from dwfconstants import *
import time
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import matlab.engine
import numpy as np

# load waveforms sdk library, expects Waveforms-SDK folder to be in current directory
current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
dwf = cdll.LoadLibrary(f"{current_directory}/Waveforms-SDK/dwf")

# device serial numbers
SCOPE1_SN = "210018B04928"
SCOPE2_SN = "210018B9D8A8"
SCOPE3_SN = "210321BB4DEB"
SUPPLY_SN = "210395B7C440"

# experiment parameters
SAMPLING_FREQUENCY = 10000 # 10 kHz sampling rate
RECORD_LENGTH = 10 # 10 second recording length
NUM_SAMPLES = RECORD_LENGTH * SAMPLING_FREQUENCY
CHANNEL_RANGE = 5.0 # +/- 5V input range
V_PLUS = 5.0 # +5V Supply
V_MINUS = -5.0 # -5V Supply
CURRENT_LIMIT = 0.4 # 10 mA to 3A range availablecs
TRIGGER_DELAY = 0.003 # 3 ms trigger delay to AD2 (empirically observed this delay during tests)
MATLAB_PROCESSING_SCRIPT = "SSA_and_ICA" # ensure matlab script is in current directory

# global digilent waveform parameters - prevent temperature drift
dwf.FDwfParamSet(DwfParamOnClose, 0) # 0 = run, 1 = stop, 2 = shutdown
dwf.FDwfParamSet(DwfParamExtFreq, 10000000) # reference clock frequency
dwf.FDwfParamSet(DwfParamFrequency, 100000000) # system clock frequency


class device:
    """
        The device class stores key identifiers for each
        digilent device (oscilloscopes and power supply).
    """
    def __init__(self, id, serial_number, description):
        self.serial_number = serial_number
        self.id = id
        self.handle = c_int()
        self.description = description
        self.acqstatus = c_int()

def open_devices():
    """
        open_devices() opens all the oscilloscopes and power
        supplies for experiments. Ensure all devices are
        powered on and connected via USB to the laptop before using
        this function.
        
        Function returns instances of the device class for the 
        three oscilloscopes and the power supply.
    """
    scope1 = None
    scope2 = None
    scope3 = None
    supply = None
    device_num = c_int(0)

    # enumerate all connected devices
    dwf.FDwfEnum(c_int(0), byref(device_num))
    print("Number of Devices Connected: " + str(device_num.value))

    # generate device objects with assigned enumeration IDs
    current_sn = create_string_buffer(16)

    for device_id in range(device_num.value):
        dwf.FDwfEnumSN(c_int(device_id), current_sn)
        sn = current_sn.value.decode('utf-8') # convert to python string for comparison
        sn = sn[3:] # remove "SN:" prefix from string
        
        if sn == SCOPE1_SN:
            scope1 = device(id=device_id, serial_number=SCOPE1_SN, description="ADP3450 Neck Radar") # neck radar scope
        elif sn == SCOPE2_SN:
            scope2 = device(id=device_id, serial_number=SCOPE2_SN, description="ADP3450 Chest Radar") # chest radar scope
        elif sn == SCOPE3_SN:
            scope3 = device(id=device_id, serial_number=SCOPE3_SN, description="AD2 ECG/PPG Sensor") # ECG/PPG scope
        elif sn == SUPPLY_SN:
            supply = device(id=device_id, serial_number=SUPPLY_SN, description="DPS3340 Power Supply") # power supply
        else:
            raise Exception(f"Unknown Device Serial Number Found: {sn}")

    # opening all devices
    for component in [scope1, scope2, scope3, supply]:
        if component == scope3:
            dwf.FDwfDeviceConfigOpen(c_int(component.id), 1, byref(component.handle))
        else:
            dwf.FDwfDeviceOpen(c_int(component.id), byref(component.handle))
        if component.handle.value == hdwfNone.value:
            szerr = create_string_buffer(512)
            dwf.FDwfGetLastErrorMsg(szerr)
            print(szerr.value)
            print("Failed to open device!")
            quit()
        print(f"{component.description} - Opened Successfully!")

    return scope1, scope2, scope3, supply

def power_supply_config(supply):
    """
        This function configures the power supply with the provided experiment parameters.
    """
    dwf.FDwfDeviceAutoConfigureSet(supply.handle, c_int(0)) # 0 = the device will only be configured when FDwf###Configure is called
    dwf.FDwfAnalogIOChannelNodeSet(supply.handle, c_int(1), c_int(0), c_double(True)) # enable output 2
    dwf.FDwfAnalogIOChannelNodeSet(supply.handle, c_int(1), c_int(1), c_double(V_MINUS)) # set v_minus
    dwf.FDwfAnalogIOChannelNodeSet(supply.handle, c_int(1), c_int(2), c_double(-1 * CURRENT_LIMIT)) 
    dwf.FDwfAnalogIOChannelNodeSet(supply.handle, c_int(2), c_int(0), c_double(True)) # enable output 3
    dwf.FDwfAnalogIOChannelNodeSet(supply.handle, c_int(2), c_int(1), c_double(V_PLUS)) # set v_plus
    dwf.FDwfAnalogIOChannelNodeSet(supply.handle, c_int(2), c_int(2), c_double(CURRENT_LIMIT))
    dwf.FDwfAnalogIOEnableSet(supply.handle, c_int(True)) # master enable on
    dwf.FDwfAnalogIOConfigure(supply.handle)
    print(f"{supply.description}: Configured!")

def adp_config(scope1, scope2):
    """
        This function configures the analog discovery pros with the provided experiment parameters.
        The analog disovery pros have a large enough sample buffer to collect all samples required
        for the experiment. Thus, the 'single' acquisition mode is used instead of 'record' for 
        the analog discovery 2.
    """
    for adp in [scope1, scope2]:
        # acquisition settings
        dwf.FDwfAnalogInChannelEnableSet(adp.handle, c_int(-1), c_int(1)) # enable all device channels
        dwf.FDwfAnalogInAcquisitionModeSet(adp.handle, acqmodeSingle)
        dwf.FDwfAnalogInChannelRangeSet(adp.handle, c_int(-1), c_double(CHANNEL_RANGE)) # set channel range for all channels
        dwf.FDwfAnalogInFrequencySet(adp.handle, c_double(SAMPLING_FREQUENCY))
        dwf.FDwfAnalogInRecordLengthSet(adp.handle, c_double(RECORD_LENGTH))
        dwf.FDwfAnalogInBufferSizeSet(adp.handle, c_int(NUM_SAMPLES))

        # trigger settings
        if adp.serial_number == SCOPE1_SN: # scope 1 used as master device
            dwf.FDwfDeviceParamSet(adp.handle, DwfParamClockMode, 1) # reference clock output on Trigger 1
            dwf.FDwfDeviceTriggerSet(adp.handle, 1, trigsrcAnalogIn) # set trigger 2 to output trigger signal
        else:
            dwf.FDwfDeviceParamSet(adp.handle, DwfParamClockMode, 2) # reference clock input on Trigger 1
            dwf.FDwfAnalogInTriggerSourceSet(adp.handle, trigsrcExternal2) # trigger signal received on trigger 2
            
        dwf.FDwfAnalogInTriggerAutoTimeoutSet(adp.handle, c_double(0)) # disable auto trigger
        dwf.FDwfAnalogInTriggerPositionSet(adp.handle, c_double(5.0))
        dwf.FDwfAnalogInConfigure(adp.handle, c_int(True), c_int(False)) # configure
        time.sleep(3)
        print(f"{adp.description}: Configured!")
        # diagnostics(adp) # uncomment if diagnostics are required to ensure device is configured correctly

def ad2_config(scope3):
    """
        This function configures the analog discovery pros with the provided experiment parameters.
        Note that the analog discovery pro collects all samples in a single buffer, however,
        the AD2 does not have a large enough buffer. Thus, the record acquisition mode is used instead of single mode.
    """
    # enable power supply
    dwf.FDwfAnalogIOChannelNodeSet(scope3.handle, c_int(0), c_int(0), c_double(True)) 
    dwf.FDwfAnalogIOChannelNodeSet(scope3.handle, c_int(0), c_int(1), c_double(3.3)) # can be changed to 5V
    dwf.FDwfAnalogIOEnableSet(scope3.handle, True) 

    # acquisition settings for analog in channels
    dwf.FDwfAnalogInChannelEnableSet(scope3.handle, c_int(-1), c_int(1))
    dwf.FDwfAnalogInAcquisitionModeSet(scope3.handle, acqmodeRecord)
    dwf.FDwfAnalogInChannelRangeSet(scope3.handle, c_int(-1), c_double(CHANNEL_RANGE)) # set channel range for all channels
    dwf.FDwfAnalogInFrequencySet(scope3.handle, c_double(SAMPLING_FREQUENCY))
    dwf.FDwfAnalogInRecordLengthSet(scope3.handle, c_double(RECORD_LENGTH))
    dwf.FDwfAnalogInBufferSizeSet(scope3.handle, c_int(16384)) # max buffer size for ad2
    
    # trigger settings
    dwf.FDwfDeviceParamSet(scope3.handle, DwfParamClockMode, 2) # reference clock input on Trigger 1
    dwf.FDwfAnalogInTriggerSourceSet(scope3.handle, trigsrcExternal2)
    dwf.FDwfAnalogInTriggerAutoTimeoutSet(scope3.handle, c_double(0)) #disable auto trigger
    dwf.FDwfAnalogInTriggerPositionSet(scope3.handle, c_double(TRIGGER_DELAY)) # trigger delay compensation - magic number for now
    dwf.FDwfAnalogInConfigure(scope3.handle, c_int(True), c_int(False)) # configure
    time.sleep(3)
    print(f"{scope3.description}: Configured!")
    # diagnostics(scope3) # uncomment if diagnostics are required to ensure device is configured correctly

def extract_scope_data(handle):
    """
        This function extracts the collected data from the analog discovery pros and returns the data
        as a python list.
    """
    data = []
    channel_num = c_int(0)
    samples_buffer = (c_double*NUM_SAMPLES)()
    dwf.FDwfAnalogInChannelCount(handle, byref(channel_num))
    for channel in range(channel_num.value):
        dwf.FDwfAnalogInStatusData(handle, channel, samples_buffer, NUM_SAMPLES)
        data.append(list(samples_buffer))
    return data

def ad2_extract(handle, ch1buffer, ch2buffer, sample_idx):
    """
        This function extracts the collected data from the analog discovery 2 using a 
        circular sample buffer. It will raise errors if samples are lost or corrupted
        during the experiment. It writes the data to the ch1 and ch2 buffers by reference
        and returns the sample index at which to start writing data again.
        
        Largely adapted from the digilent waveforms analogin_record example scripts.
    """
    available = c_int()
    lost = c_int()
    corrupted = c_int()
    
    dwf.FDwfAnalogInStatusRecord(handle, byref(available), byref(lost), byref(corrupted))
    sample_idx += lost.value
    sample_idx %= NUM_SAMPLES

    if lost.value :
        raise Exception("AD2 samples were lost during reocrding")
    if corrupted.value :
        raise Exception("AD2 samples were corrupted during recording")

    iBuffer = 0
    while available.value>0:
        cSamples = available.value
        # using circular sample buffer
        if sample_idx + available.value > NUM_SAMPLES:
            cSamples = NUM_SAMPLES - sample_idx
        dwf.FDwfAnalogInStatusData2(handle, c_int(0), byref(ch1buffer, sizeof(c_double)*sample_idx), c_int(iBuffer), c_int(cSamples)) # get channel 1 data
        dwf.FDwfAnalogInStatusData2(handle, c_int(1), byref(ch2buffer, sizeof(c_double)*sample_idx), c_int(iBuffer), c_int(cSamples)) # get channel 2 data
        iBuffer += cSamples
        available.value -= cSamples
        sample_idx += cSamples
        sample_idx %= NUM_SAMPLES
    
    return sample_idx

def run_experiment(scope1, scope2, scope3):
    """
        This function handles running the experiment and returning
        all the data collected as a pandas dataframe.
        
        To check the acquisition states are as expected, you can use the following lines.
        The acqstatus values and what they mean are documented in the waveforms sdk manual.
        dwf.FDwfAnalogInStatus(scope1.handle, 1, byref(scope1.acqstatus))
        print(scope1.acqstatus.value)
    """
    
    print("Acquisition Starting...")
    # initialize data recording buffers for analog discovery 2
    ch1buffer = (c_double*NUM_SAMPLES)()
    ch2buffer = (c_double*NUM_SAMPLES)()
    sample_idx = 0 # index for keeping track of AD2 recording process
    
    # configure scopes to start reading data upon trigger signal (master configured last)
    dwf.FDwfAnalogInConfigure(scope2.handle, 1, 1)
    dwf.FDwfAnalogInConfigure(scope3.handle, 1, 1)
    time.sleep(2)
    dwf.FDwfAnalogInConfigure(scope1.handle, 1, 1) # master configured last - sends a trigger signal to all devices
        
    while True:
        sample_idx = ad2_extract(scope3.handle, ch1buffer, ch2buffer, sample_idx)
        
        dwf.FDwfAnalogInStatus(scope1.handle, 1, byref(scope1.acqstatus))
        dwf.FDwfAnalogInStatus(scope2.handle, 1, byref(scope2.acqstatus))
        dwf.FDwfAnalogInStatus(scope3.handle, 1, byref(scope3.acqstatus))
        
        if scope1.acqstatus.value == DwfStateDone.value and scope2.acqstatus.value == DwfStateDone.value and scope3.acqstatus.value == DwfStateDone.value:
            print("Acquisition done!")
            break
    
        # script downtime for ad2 sample buffer to fill up. 5000 samples is a magic number.
        # if ad2 samples are getting lost or corrupted and causing errors, lower the sleep time
        # by reducing the number of samples collected each time.
        time.sleep(5000.0 / SAMPLING_FREQUENCY)
    
    neck_radar_data = extract_scope_data(scope1.handle)
    chest_radar_data = extract_scope_data(scope2.handle)
    
    # adjust ad2 channel data circular sample buffers according to sample index
    # for correctly ordered time series data
    if sample_idx != 0 :
        ch1buffer = ch1buffer[sample_idx:]+ch1buffer[:sample_idx]
        ch2buffer = ch2buffer[sample_idx:]+ch2buffer[:sample_idx]
    
    # post-processing matlab scripts expect columns in this order
    experiment_data = {}
    experiment_data["ECG"] = list(ch2buffer)
    experiment_data["PPG"] = list(ch1buffer)
    experiment_data["Left Neck - I"] = neck_radar_data[1]
    experiment_data["Left Neck - Q"] = neck_radar_data[0]
    experiment_data["Right Neck - I"] = neck_radar_data[3]
    experiment_data["Right Neck - Q"] = neck_radar_data[2]

    experiment_data["Left Chest - I"] = chest_radar_data[1]
    experiment_data["Left Chest - Q"] = chest_radar_data[0]
    experiment_data["Right Chest - I"] = chest_radar_data[3]
    experiment_data["Right Chest - Q"] = chest_radar_data[2]
    
    final_data = pd.DataFrame(experiment_data)
    
    return final_data

def plot_results(experiment_data, save=False, filepath=None):
    """
        This function handles all plotting the raw data
        and radar complex signal IQ plots. If 'save' is set
        to True, the plots will be saved to the given filepath.
    """
    
    fig, axs = plt.subplots(3, 1, figsize=(18, 8))
    
    # Plot Neck data
    axs[0].plot(experiment_data["Left Neck - I"], label="Left Neck - I")
    axs[0].plot(experiment_data["Left Neck - Q"], label="Left Neck - Q")
    axs[0].plot(experiment_data["Right Neck - I"], label="Right Neck - I")
    axs[0].plot(experiment_data["Right Neck - Q"], label="Right Neck - Q")
    axs[0].set_title("Neck Signals")
    axs[0].set_xlabel("Sample Number")
    axs[0].legend(loc='upper right')

    # Plot Chest data
    axs[1].plot(experiment_data["Left Chest - I"], label="Left Chest - I")
    axs[1].plot(experiment_data["Left Chest - Q"], label="Left Chest - Q")
    axs[1].plot(experiment_data["Right Chest - I"], label="Right Chest - I")
    axs[1].plot(experiment_data["Right Chest - Q"], label="Right Chest - Q")
    axs[1].set_title("Chest Signals")
    axs[1].set_xlabel("Sample Number")
    axs[1].legend(loc='upper right')

    # # Plot ECG/PPG data
    axs[2].plot(experiment_data["ECG"], label="ECG")
    axs[2].plot(experiment_data["PPG"], label="PPG")
    axs[2].set_title("ECG/PPG Signals")
    axs[2].set_xlabel("Sample Number")
    axs[2].legend(loc='upper right')
    
    # Adjust layout and display
    plt.tight_layout()
    if save:
        plt.savefig(f"{filepath}/Raw_data.png")
        plt.close()
    else:
        plt.show()

    fig, axs = plt.subplots(2, 2, figsize=(9, 9))
    step_size = 100 # avoid overflow errors when plotting large amounts of data
    
    # Left Neck Radar IQ Plot
    axs[0, 0].plot(experiment_data["Left Neck - I"][::step_size], experiment_data["Left Neck - Q"][::step_size], label="Left Neck - IQ")
    axs[0, 0].set_title("Left Neck Radar IQ Plot")
    axs[0, 0].set_xlabel("I")
    axs[0, 0].set_ylabel("Q")
    axs[0, 0].legend()
    
    # Right Neck Radar IQ Plot
    axs[0, 1].plot(experiment_data["Right Neck - I"][::step_size], experiment_data["Right Neck - Q"][::step_size], label="Right Neck - IQ")
    axs[0, 1].set_title("Right Neck Radar IQ Plot")
    axs[0, 1].set_xlabel("I")
    axs[0, 1].set_ylabel("Q")
    axs[0, 1].legend()

    # Left Chest Radar IQ Plot
    axs[1, 0].plot(experiment_data["Left Chest - I"][::step_size], experiment_data["Left Chest - Q"][::step_size], label="Left Chest - IQ")
    axs[1, 0].set_title("Left Chest Radar IQ Plot")
    axs[1, 0].set_xlabel("I")
    axs[1, 0].set_ylabel("Q")
    axs[1, 0].legend()
    
    # Right Chest Radar IQ Plot
    axs[1, 1].plot(experiment_data["Right Chest - I"][::step_size], experiment_data["Right Chest - Q"][::step_size], label="Right Chest - IQ")
    axs[1, 1].set_title("Right Chest Radar IQ Plot")
    axs[1, 1].set_xlabel("I")
    axs[1, 1].set_ylabel("Q")
    axs[1, 1].legend()
    
    # Adjust layout and display
    plt.tight_layout()
    if save:
        plt.savefig(f"{filepath}/Raw_IQ_Plots.png")
        plt.close()
    else:
        plt.show()

def diagnostics(scope):
    """
        Helper function to check various device parameters and ensure
        commands are being executed as expected. Use for debugging code
        as waveforms sdk functions tend to not work as expected all the time.
    """
    
    max_buffer = c_int()
    min_buffer = c_int()
    dwf.FDwfAnalogInBufferSizeInfo(scope.handle, byref(min_buffer), byref(max_buffer))
    print(f"Max buffer size ({scope.description}): {max_buffer.value}")
    print(f"Min buffer size ({scope.description}): {min_buffer.value}")
    
    current_buffer_size = c_int()
    dwf.FDwfAnalogInBufferSizeGet(scope.handle, byref(current_buffer_size))
    print(f"Current buffer size ({scope.description}): {current_buffer_size.value}")

    trigger_position = c_double()
    dwf.FDwfAnalogInTriggerPositionGet(scope.handle, byref(trigger_position))
    print(f"Current trigger position ({scope.description}): {trigger_position.value}")
    
    auto_configure = c_int()
    dwf.FDwfDeviceAutoConfigureGet(scope.handle, byref(auto_configure))
    print(f"Current autoconfigure setting ({scope.description}): {auto_configure.value}")
    
    device_trigger = c_int()
    dwf.FDwfDeviceTriggerGet(scope.handle, c_int(0), byref(device_trigger))
    print(f"Current device trigger setting (T1) ({scope.description}): {device_trigger.value}")
    
    dwf.FDwfDeviceTriggerGet(scope.handle, c_int(1), byref(device_trigger))
    print(f"Current device trigger setting (T2) ({scope.description}): {device_trigger.value}")
    
    analogin_trigger = c_int()
    dwf.FDwfAnalogInTriggerSourceGet(scope.handle, byref(analogin_trigger))
    print(f"Current AnalogIn trigger setting ({scope.description}): {analogin_trigger.value}")
    
    status = c_int()
    dwf.FDwfAnalogInStatus(scope.handle, c_int(1), byref(status))
    print(f"Current scope status: ({scope.description}): {device_trigger.value}")
    
    acqmode = c_int()
    dwf.FDwfAnalogInAcquisitionModeGet(scope.handle, byref(acqmode))
    print(f"Current acquisition mode: ({scope.description}): {acqmode.value}")
    
    offset = c_double()
    dwf.FDwfAnalogInChannelOffsetGet(scope.handle, c_int(3), byref(offset))
    print(f"Channel offset (Channel 4 as an example): ({scope.description}): {offset.value}")
    
def matlab_post_processing(filepath):
    """
        This function starts a matlab engine and executes post-processing
        on saved experiment data. Post-processing pipeline includes
        Singular Spectrum Analysis and Independent Components Analysis
        in the SSA_and_ICA.m file.
        
        The matlab file that will be run is specified in the eng.eval line.
        Ensure the file is in the python script's directory.
    """
    eng = matlab.engine.start_matlab() # start matlab engine
    eng.cd(os.path.dirname(os.path.abspath(__file__)), nargout=0)
    
    eng.workspace['data_filepath'] = filepath
    eng.workspace['ECG'] = 2
    eng.workspace['PPG'] = 3

    # post-processing neck data
    print("Matlab neck radar data processing...")
    eng.workspace['I1'] = 4
    eng.workspace['Q1'] = 5
    eng.workspace['Radar1'] = 'Left_Neck'
    eng.workspace['I2'] = 6
    eng.workspace['Q2'] = 7
    eng.workspace['Radar2'] = 'Right_Neck'
    eng.workspace['position'] = 'Neck'
    eng.eval(MATLAB_PROCESSING_SCRIPT, nargout=0)

    # post-processing chest data
    print("Matlab chest radar processing...")
    eng.workspace['I1'] = 8
    eng.workspace['Q1'] = 9
    eng.workspace['Radar1'] = 'Left_Chest'
    eng.workspace['I2'] = 10
    eng.workspace['Q2'] = 11
    eng.workspace['Radar2'] = 'Right_Chest'
    eng.workspace['position'] = 'Chest'
    eng.eval(MATLAB_PROCESSING_SCRIPT, nargout=0)
    
    eng.quit() # close matlab engine

    
def main():
    scope1, scope2, scope3, supply = open_devices()
    power_supply_config(supply)
    adp_config(scope1, scope2)
    ad2_config(scope3)

    # Menu to run several experiments
    while True:
        print("\nMenu:")
        print("1. Run Experiment")
        print("2. Quit")
        choice = input("Enter your choice: ")
        
        while choice != '1' and choice != '2':
            print("Invalid Choice, Try Again...")
            choice = input("Enter your choice: ")
        
        if choice == '1':
            experiment_data = run_experiment(scope1, scope2, scope3)
            plot_results(experiment_data)
            save = input("Save experiment data and run matlab post-processing? (y/n)?: ")
            while save != 'y' and save != 'n':
                print("Invalid Choice, Try Again...")
                save = input("Save experiment data and run matlab post-processing? (y/n)?: ")
            if save == 'y':
                filename = input("Enter filename for data: ")
                os.makedirs(f"{current_directory}/Data/{filename}", exist_ok=True)
                filepath = f"{current_directory}/Data/{filename}/{filename}.csv"
                experiment_data.to_csv(filepath, index=True)
                plot_results(experiment_data, save=True, filepath=f"{current_directory}/Data/{filename}")
                print("Data/Figs saved in 'Data' folder")
                print("Starting Matlab Post-Processing...")
                matlab_post_processing(filepath)
                print("Save Complete!")
            elif save == 'n':
                print("Data not saved")
        elif choice == '2':
            dwf.FDwfDeviceCloseAll() # close all devices
            break
    
    print("Session Terminated!")

if __name__ == '__main__':
    main()