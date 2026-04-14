"""
Functions to control the powermeters
"""
import pyvisa as visa
from ThorlabsPM100 import ThorlabsPM100

class ThorlabsPowermeter:
    
    def __init__(self, address, wavelength=1550):
        rm = visa.ResourceManager()
        #devs = rm.list_resources()
        powerMeter = rm.open_resource(address)
        powerMeter.write('CONF:POW')
        
        self.pm = ThorlabsPM100(inst=powerMeter)
        self.set_correction_wavelength(wavelength)
        
    def set_correction_wavelength(self, wavelength):
        self.pm.sense.correction.wavelength = wavelength
        
    def read(self):
        return self.pm.read
        