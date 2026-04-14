import pyvisa as visa
from .qontrol import Qontrol
from .powermeter import ThorlabsPowermeter
# from .qutag import QuToolsTimeTagger, record_counts, record_coincidences
#from .swabian import SwabianTimeTagger
from .laser import CobriteLaser, SantecLaser

def list_devices():
    rm = visa.ResourceManager()
    devs = rm.list_resources()
    print(devs)
