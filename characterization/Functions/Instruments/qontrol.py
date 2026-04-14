"""
Functions to control the Qontrol modules
"""
import qontrol

class Qontrol:
    
    def __init__(self, COM=4):
        self.q = qontrol.QXOutput(serial_port_name="COM" + str(COM), response_timeout=1)
        
    def reset(self):
        self.q.v[:] = 0
    
    def read_i(self, channel=None):
        if channel is not None:
            return self.q.i[channel]
        else:
            return self.q.i
    
    def read_v(self, channel=None):
        if channel is not None:
            return self.q.v[channel]
        else:
            return self.q.v
        
    def set_maxv(self, channel, max_voltage):
        self.q.vmax[channel] = max_voltage
        print("Max voltage for channel {} is {} V".format(channel, self.q.vmax[channel]))

    def set_maxi(self, channel, max_current):
        self.q.imax[channel] = max_current
        print("Max current for channel {} is {} A".format(channel, self.q.imax[channel]))
    
    def set_v(self, channel, voltage):
        # if voltage < 0 or voltage > 6:
        #     raise ValueError("Voltage must be between 0 and 6 V")
        self.q.v[channel] = voltage
        
    def close(self):
        self.q.close()
    


        