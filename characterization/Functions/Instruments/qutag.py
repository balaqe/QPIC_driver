# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 14:58:03 2025

@author: FTNK-LocalAdm
"""
#%% Imports 
from Functions.BasicFunctions import *
from QtagCodes.TimeTaggerCallback.QuTAG_MC import QuTAG
import sys

#%% Connect time tagger

class QuToolsTimeTagger():

    def __init__(self, expTime=1000, coincWin=200, channel_delay=None, delay_ps=0):
        # Create a QuTag object and initialize it

        folder_path = r"\\ait-pdfs.win.dtu.dk\Services\ELEC\comm\QuantumLab\chip_exps\Lab_Codes\QtagCodes\TimeTaggerCallback"
        sys.path.append(folder_path)
        os.chdir(folder_path)
        from QuTAG_MC import QuTAG  

        self.qutag = QuTAG(-1)
        self.qutag.init()

        # Get the version and Timebase
        _deviceType = self.qutag.getVersion()
        print("Device Type: ", _deviceType)
        _timebase = self.qutag.getTimebase()
        print("Time base: ", _timebase)

        # Update the number of available channels
        self.qutag.TDC_QUTAG_CHANNELS = self.qutag.qutools_dll.TDC_getChannelCount()

        # Set the buffer size after initialization
        self.qutag.setTimestampBufferSize(self.qutag._bufferSize)

        # Set the parameters
        self.expTime = expTime
        self.coincWin = coincWin

        self.set_exposure(self.expTime)
        self.set_coincidence_window(self.coincWin)
        self.set_delay(channel_delay, delay_ps)
        
        self.qutag.freezeBuffers(False)
        self.qutag.enableMarkers([4])

        _, self.coincWin, self.expTime = self.qutag.getDeviceParams()
        print("Coincidence window:", self.coincWin, "bins, Exposure time: ", self.expTime, "ms")

        #" Making sure Start is OFF "
        #self.qutag.enableChannels(0, int('11111111', 2))

    def set_exposure(self, exp_time_ms):
        self.expTime = exp_time_ms
        self.qutag.setExposureTime(self.expTime)

    def set_coincidence_window(self, coinc_win):
        self.coincWin = coinc_win
        self.qutag.setCoincidenceWindow(self.coincWin)

    def set_delay(self, channel_delay, delay_ps):
        if channel_delay is not None:
            if isinstance(channel_delay, int):
                self.qutag.setChannelDelay(channel_delay, delay_ps)
            else:
                for i in range(len(channel_delay)):
                    self.qutag.setChannelDelay(channel_delay[i], delay_ps[i])

    def get_counts(self, channel=None):
        self.qutag.freezeBuffers(True)
        ret, data, update = self.qutag.getCoincCounters()
        if channel is not None:
            return data[channel]
        else:
            return data

def record_counts(qutag:QuTAG, channels, num_measurements, folder=None, filename=None):
    if folder is None:
        folder = select_folder()
    print("Estimated acquisition time: ", num_measurements*qutag.expTime/1000, "s")
    counts_arr = [[] for _ in range(len(channels))]
    for i in range(num_measurements):
        for i in range(len(channels)):
            counts_arr[i].append(qutag.get_counts(channels[i]))
        time.sleep(qutag.expTime/1000)  
    
    # concatenate all the elements of counts_arr in a single array¨
    channels_arr = np.concatenate([channels[i]*np.ones(len(counts_arr[i])) for i in range(len(channels))])
    counts_arr = np.concatenate(counts_arr)


    if filename is None: filename = 'Coincidences'
    filepath = folder + "\\" + filename

    i = 0
    while os.path.exists(filepath + '.txt'):
        filepath = folder + "\\" + filename + "_" + str(i)
        i += 1
    np.savetxt(filepath+'.txt', np.column_stack((channels_arr, counts_arr)), delimiter=' ; ', header='Channel ; Coincidences')
    # print("Saved data to: ", filepath)


        
# Record the coincidence counts, channel 38 = detectors 3/4
def record_coincidences(qt:QuTAG, num_measurements, channel=38, folder=None, filename=None):
    if folder is None: 
        folder = select_folder()
    # print("Estimated acquisition time: ", num_measurements*qutag.expTime/1000, "s")
    coincidences_arr = []
    time.sleep(1)
    for i in range(num_measurements):
        coincidences_arr.append(qt.get_counts(channel))
        time.sleep(qt.expTime/1000) 
 
    # concatenate all the elements of counts_arr in a single array
    channel_arr = channel * np.ones(len(coincidences_arr))

    if filename is None: filename = 'Coincidences'
    filepath = folder + "\\" + filename

    i = 0
    while os.path.exists(filepath + '.txt'):
        filepath = folder + "\\" + filename + "_" + str(i)
        i += 1
    np.savetxt(filepath+'.txt', np.column_stack((channel_arr, coincidences_arr)), delimiter=' ; ', header='Channel ; Coincidences')
    # print("Saved data to: ", filepath)

    return coincidences_arr


# def test_counts(qt:QuTAG):
#         qt.freezeBuffers(True)
#         ret, data, update = qt.getCoincCounters()
#         print(data)


#%%
if __name__ == '__main__':
    qt = QuToolsTimeTagger(expTime=1000, coincWin=200, channel_delay=2, delay_ps=31840)
    #print(qt.get_counts()) 
    # record_counts(qt, [3,4], 10*1)
    record_coincidences(qt, num_measurements=15, channel=35)


# %%
