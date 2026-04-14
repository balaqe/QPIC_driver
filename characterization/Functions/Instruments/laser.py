"""
Functions to control the laser
"""
import serial
import time
import pyvisa as visa

class PhotoneticsLaser:
    
    def __init__(self):
        pass
    

class SantecLaser:
    
    def __init__(self, GPIB=28, timeout=3000):
        self.laser = None
        self.gpib = GPIB
        #self.timeout = timeout
        self.power_dbm = None
        self.wavelength_nm = None

    def connect(self):
        self.laser =visa.ResourceManager().open_resource('GPIB0::' + str(self.gpib) + '::INSTR')

    def turn_on(self):
        self.laser.write('LO')

    # Round the wavelength to the nearest 0.01 nm, otherwise the laser can get stuck
    def set_wavelength(self, wl):
        self.laser.write(f'WA{round(wl, 2)}')

    def set_power(self, power):
        self.laser.write(f'OP{power}')  

    def turn_off(self):
        self.laser.write('LF')

    def close(self):
        self.laser.close()
    

class CobriteLaser:

    def __init__(self, COM:int, baudrate=115200):
        self.com_port = "COM" + str(COM)
        self.baudrate = baudrate
        self.ser = None
        self.sep = ";" 

        self.power_dbm = None
        self.wavelength_nm = None 

    def connect(self):
        try:
            self.ser = serial.Serial(self.com_port, self.baudrate, timeout=2)
            self.ser.reset_input_buffer()
            self.ser.reset_output_buffer()
            print(f"Connected to COM {self.com_port}")
            self.query("pass IDP")
            idn = self.query("*idn?")
            print(f"Device ID: {idn}")
        except Exception as e:
            print(f"Connection Error: {e}")

    def query(self, command):
        if not self.ser: return
        full_cmd = f"{command}{self.sep}"
        self.ser.write(full_cmd.encode('utf-8'))
        reply = ""
        start_time = time.time()
        while self.sep not in reply:
            if (time.time() - start_time) > 5: 
                break
            if self.ser.in_waiting > 0:
                reply += self.ser.read(self.ser.in_waiting).decode('utf-8')
        clean_reply = reply.replace(self.sep, "").strip()
        #print(f"TX: {command} -> RX: {clean_reply}")
        return clean_reply

    def set_wavelength(self, wavelength_nm, card=1, port=1):
        self.wavelength_nm = wavelength_nm
        self.query(f"wav {card},1,{port},{wavelength_nm}")

        # print("Tuning... please wait.", end='\r')
        # self.query(f"bwai {card},1,{port}")
        # print("Laser is stable and ON.          ")

    def set_power(self, power_dbm, card=1, port=1):
        self.power_dbm = power_dbm
        self.query(f"pow {card},1,{port},{power_dbm}")

        # print("Tuning... please wait.", end='\r')
        # self.query(f"bwai {card},1,{port}")
        # print("Laser is stable and ON.          ")

    def turn_on(self, card=1, port=1):
        self.query(f"stat {card},1,{port},1")

    def turn_off(self, card=1, port=1):
        self.query(f"stat {card},1,{port},0")

    def close(self):
        if self.ser:
            self.ser.close()
            print("Connection closed.")