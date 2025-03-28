"""
   DWF Python Example
   Author:  Digilent, Inc.
   Revision:  2023-12-08

   Requires:                       
       Python 2.7, 3
"""

from ctypes import *
from dwfconstants import *
import time
import sys

if sys.platform.startswith("win"):
    dwf = cdll.dwf
elif sys.platform.startswith("darwin"):
    dwf = cdll.LoadLibrary("/Library/Frameworks/dwf.framework/dwf")
else:
    dwf = cdll.LoadLibrary("libdwf.so")

hdwf = c_int()
vpp = c_double()
vpn = c_double()
pcbTemperature = c_double()
dieTemperature = c_double()
usbVoltage = c_double()
usbCurrent = c_double()
auxVoltage = c_double()
auxCurrent = c_double()
suppliesLimit = c_double()
dsts = c_double()

version = create_string_buffer(16)
dwf.FDwfGetVersion(version)
print("DWF Version: "+str(version.value))
szerr = create_string_buffer(512)

#open device
print("Opening first device")

if dwf.FDwfDeviceOpen(c_int(-1), byref(hdwf)) == 0 or hdwf.value == hdwfNone.value:
    dwf.FDwfGetLastErrorMsg(szerr)
    print("Failed to open device:")
    print(szerr.value.decode())
    quit()

dwf.FDwfGetLastErrorMsg(szerr)
if len(szerr.value) != 0:
    print("Warning:")
    print(szerr.value.decode())


print("Device temperature and USB/AUX supply voltage and current")

# 10 times, once per second
for i in range(1, 11):
    # wait between readings
    time.sleep(1)
    # fetch analog IO status from device
    if dwf.FDwfAnalogIOStatus(hdwf) == 0 :
        dwf.FDwfGetLastErrorMsg(szerr)
        print("Communication error:")
        print(szerr.value.decode())
        break;
    # get system monitor readings
    dwf.FDwfAnalogIOChannelNodeStatus(hdwf, c_int(0), c_int(1), byref(vpp))
    dwf.FDwfAnalogIOChannelNodeStatus(hdwf, c_int(1), c_int(1), byref(vpn))
    dwf.FDwfAnalogIOChannelNodeStatus(hdwf, c_int(2), c_int(0), byref(pcbTemperature))
    dwf.FDwfAnalogIOChannelNodeStatus(hdwf, c_int(2), c_int(1), byref(dieTemperature))
    dwf.FDwfAnalogIOChannelNodeStatus(hdwf, c_int(2), c_int(2), byref(usbVoltage))
    dwf.FDwfAnalogIOChannelNodeStatus(hdwf, c_int(2), c_int(3), byref(usbCurrent))
    dwf.FDwfAnalogIOChannelNodeStatus(hdwf, c_int(2), c_int(4), byref(auxVoltage))
    dwf.FDwfAnalogIOChannelNodeStatus(hdwf, c_int(2), c_int(5), byref(auxCurrent))
    dwf.FDwfAnalogIOChannelNodeStatus(hdwf, c_int(3), c_int(0), byref(suppliesLimit))
    dwf.FDwfAnalogIOChannelNodeStatus(hdwf, c_int(15), c_int(0), byref(dsts))
    
    if dsts.value == 1: print("Device is ON")
    else: print("Device is OFF")
    print("Supplies Limit: " + str(suppliesLimit.value) + " W")
    print("Positive Supply: " + str(round(vpp.value,3)) + " V")
    print("Negative Supply: " + str(round(vpn.value,3)) + " V")
    print("PCB Temperature: " + str(round(pcbTemperature.value,2)) + "*C")
    print("Die Temperature: " + str(round(dieTemperature.value,2)) + "*C")
    print("USB:\t" + str(round(usbVoltage.value,3)) + "V\t" + str(round(usbCurrent.value,3)) + "A")
    print("AUX:\t" + str(round(auxVoltage.value,3)) + "V\t" + str(round(auxCurrent.value,3)) + "A")
    print("")
    
#close the device
dwf.FDwfDeviceClose(hdwf)
