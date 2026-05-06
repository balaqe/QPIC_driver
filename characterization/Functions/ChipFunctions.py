import numpy as np
import cmath
import time
import os
import json
from scipy.optimize import curve_fit
from .BasicFunctions import *
import matplotlib.pyplot as plt

"A: amplitude"
"B: amplitude"
"omega: frequency"
"phi_0: phase offset"
"a, b, c: from the I-V fit"

"P: -phi_0 -phi"
"rho_0: omega * c"
"rho_1: omega * b"
"rho_2: omega * a"

"phi: phase you want to set"


class Chip():
    def __init__(self, qontrol, folder=None):
        print(f"Chip initialized, mzi_dict set: {id(self)}")
        self.folder = folder or select_folder()
        self.qontrol = qontrol
        self.mzi_dict = {}
        self.ext_phase_dict = {}
        self.load_params(folder)
    
    def set_config(self, data):
        for m_name, data_point in data.items():
            print(f"type of datapoint: {type(data_point)}")
            if 'shifter1' in data_point: # Check if it's an MZI
                ps = self.mzi_dict[m_name].shifter[0]
                ps.set_phase(data_point['shifter1']*np.pi)
            if 'shifter2' in data_point:
                ps = self.mzi_dict[m_name].shifter[1]
                ps.set_phase(data_point['shifter2']*np.pi)

            if 'shifter1' not in data_point and 'shifter2' not in data_point:
                ps = self.ext_phase_dict[m_name]
                ps.set_phase(data_point*np.pi)

        # for key, mzi in phases.items():

        #     ps = PhaseShifter(self, mzi)
        #     if not os.path.exists(ps.folder + 'data.json'):
        #         raise FileNotFoundError(f'data.json file not found in {ps.folder}')
        #     else:
        #         ps.set_phase(phases[mzi])

    def load_config(self, file_path, keys_arr=None):

        def get_nested(data, keys):
            if not keys or keys is None:
                return data
            else:
                return get_nested(data[keys[0]], keys[1:])
        #
        # check if keys is a list
        if keys_arr is not None and not isinstance(keys_arr, list):
            raise TypeError("keys_arr must be a list")

        with open(file_path, "r") as f:
            datafile = json.load(f)

        chip_config = get_nested(datafile, keys_arr)    
        chip_config = {key: value * np.pi for key, value in chip_config.items()}
        self.set_config(chip_config)    

    def add_mzi(self, mzi: MZI):
        self.mzi_dict[mzi.name] = mzi

    def add_ext_phase(self, phase_shifter: PhaseShifter):
        self.ext_phase_dict[phase_shifter.name] = phase_shifter

    def get_dictionary(self):
        res = {}
        for _, mzi in self.mzi_dict.items():
            res[mzi.name] = {
                "shifter1": mzi.shifter[0].get_master_params(),
                "shifter2": mzi.shifter[1].get_master_params()
            }
            mzi.shifter[0].save_parameters()
            mzi.shifter[1].save_parameters()

        for _, shifter in self.ext_phase_dict.items():
            res[shifter.name] = shifter.get_master_params()
            shifter.save_parameters()

        return res

    def save_parameters(self):
        save_path = os.path.join(self.folder, "data_master.json")
        save_to_json(save_path, name="clements", data=self.get_dictionary())

    def load_params(self, file_path):
        if os.path.exists(file_path):
            if "data_master.json" not in file_path:
                file_path = os.path.join(file_path, "data_master.json")
            with open(file_path, "r") as f:
                data = json.load(f)
        
            self.mzi_dict = {}

            for m_name, data_point in data.items():
                print(f"m_name found: {m_name}")
                if "shifter1" in data_point: # Check if it's an MZI
                    shifter1 = PhaseShifter(chip=self, name=data_point["shifter1"]["name"])
                    shifter2 = PhaseShifter(chip=self, name=data_point["shifter2"]["name"])

                    shifter1.set_params(data_point["shifter1"])
                    shifter2.set_params(data_point["shifter2"])

                    mzi = MZI(chip=self, name=m_name, shifter1=shifter1, shifter2=shifter2)
                    self.mzi_dict[mzi.name] = mzi
                else:
                    self.ext_phase_dict[m_name] = PhaseShifter(chip=self, name=m_name)

class MZI():
    def __init__(self, chip, name, shifter1: PhaseShifter, shifter2: PhaseShifter, folder=None):
        self.name = name
        self.chip = chip
        self.folder = folder or self.chip.folder + f'/MZI_{name}/'
        create_folder(self.folder)
        self.shifter = [shifter1, shifter2]


class PhaseShifter():
    
    def __init__(self, chip, name, folder=None):
        self.chip = chip
        self.folder = folder or self.chip.folder + f'/MZI_{name}/'
        create_folder(self.folder)
        
        self.name = name
        
        self.in_port = None
        self.out_port = None
        self.channel = None

        self.phase_flip = False

        self.temperature = None

        self.chip_config = None

        self.temperature = None
        
        self.volt_0 = None
        self.volt_pi = None
        self.volt_pi2 = None
        self.volt_3pi2 = None
        
        self.visibility = None
        
        self.rho0 = None
        self.rho1 = None
        self.rho2 = None
        
        self.A = None
        self.B = None
        self.omega = None
        self.phi_0 = None
        
        self.a = None
        self.b = None
        self.c = None
        
        self.voltage_arr = None
        self.current_arr = None
        self.opt_voltage_arr = None
        self.opt_power_arr = None
        
        # Initialization of the attributes if one file already exists
        # self.load_parameters()

    def set_params(self, input):
        self.in_port = input["in_port"]
        self.out_port = input["out_port"]
        self.channel = input["channel"]

        self.phase_flip = input["phase_flip"]

        self.temperature = input["temperature"]

        self.chip_config = input["chip_config"]
        
        self.volt_0 = input["volt_0"]
        self.volt_pi = input["volt_pi"]
        self.volt_pi2 = input["volt_pi2"]
        self.volt_3pi2 = input["volt_3pi2"]
        
        self.visibility = input["visibility"]
        
        self.rho0 = input["rho0"]
        self.rho1 = input["rho1"]
        self.rho2 = input["rho2"]
        
        self.A = input["A"]
        self.B = input["B"]
        self.omega = input["omega"]
        self.phi_0 = input["phi_0"]
        
        self.a = input["a"] 
        self.b = input["b"] 
        self.c = input["c"] 
        
    def get_params_dict(self):
            return {
                "name": self.name,
                "in_port": self.in_port,
                "out_port": self.out_port,
                "channel": self.channel,
                "phase_flip": self.phase_flip,
                "temperature": self.temperature,
                "chip_config": self.chip_config,
                "volt_0": self.volt_0,
                "volt_pi": self.volt_pi,
                "volt_pi2": self.volt_pi2,
                "volt_3pi2": self.volt_3pi2,
                "visibility": self.visibility,
                "rho0": self.rho0,
                "rho1": self.rho1,
                "rho2": self.rho2,
                "A": self.A,
                "B": self.B,
                "omega": self.omega,
                "phi_0": self.phi_0,
                "a": self.a,
                "b": self.b,
                "c": self.c,
                "voltage_arr": self.voltage_arr,
                "current_arr": self.current_arr,
                "opt_voltage_arr": self.opt_voltage_arr,
                "opt_power_arr": self.opt_power_arr
            }
    
    def get_master_params(self):
            return {
                "name": self.name,
                "in_port": self.in_port,
                "out_port": self.out_port,
                "channel": self.channel,
                "chip_config": self.chip_config,
                "temperature": self.temperature,
                "phase_flip": self.phase_flip,
                "volt_0": self.volt_0,
                "volt_pi": self.volt_pi,
                "volt_pi2": self.volt_pi2,
                "volt_3pi2": self.volt_3pi2,
                "visibility": self.visibility,
                "rho0": self.rho0,
                "rho1": self.rho1,
                "rho2": self.rho2,
                "A": self.A,
                "B": self.B,
                "omega": self.omega,
                "phi_0": self.phi_0,
                "a": self.a,
                "b": self.b,
                "c": self.c
            }
            
    def load_parameters(self):
        file_path = self.folder + 'data.json'
        if os.path.exists(file_path):
            with open(file_path, "r") as f:
                datafile = json.load(f)
            for key, value in datafile.items():
                if hasattr(self, key): 
                    setattr(self, key, value)

    def save_parameters(self):
        for key, value in self.get_params_dict().items():
            save_to_json(self.folder + 'data.json', key, value)
        
    def run_electrical_sweep(self, start_v=0, stop_v=6.1, end_v=0.1, delay=0.02):
        measured_voltage_arr = []
        measured_current_arr = []
        for voltage in np.arange(start_v, stop_v, end_v):
            print('%.2f' % (voltage), end='\r')
            self.chip.qontrol.set_v(self.channel, float(voltage))
            time.sleep(delay)
            measured_voltage_arr.append(self.chip.qontrol.read_v(self.channel))
            measured_current_arr.append(self.chip.qontrol.read_i(self.channel))
        # set all channels to zero
        # qontrol.reset()
        # qontrol.close()
        self.voltage_arr = np.asarray(measured_voltage_arr)
        self.current_arr = np.asarray(measured_current_arr)
        
    def run_electrical_fit(self):
        self.a, self.b, self.c = np.polyfit(self.voltage_arr, self.current_arr, 2)
    
    def plot_electrical_fit(self):
        plt.plot(self.voltage_arr, self.a * self.voltage_arr**2 + self.b * self.voltage_arr + self.c, c='b', label = 'Fit', zorder=-1)
        plt.scatter(self.voltage_arr, self.current_arr, s=15, facecolors='w', edgecolors='orchid', label = 'Data') 
        plt.xlabel('V(V)')
        plt.ylabel('i(mA)')
        plt.title('MZI '+ self.name)
        plt.grid(linestyle=':', alpha=0.5)
        plt.savefig(self.folder+'MZI_'+self.name+'_electrical.png', dpi=300)
        plt.show()
        plt.close()
        
    def plot_electrical_power(self):
        i_fit = self.a * self.voltage_arr**2 + self.b * self.voltage_arr + self.c
        power = i_fit * self.voltage_arr
        plt.scatter(self.voltage_arr, power)
        plt.xlabel('V(V)')
        plt.ylabel('P(mW)')
        plt.title('MZI '+ self.name)
        plt.grid()
        plt.close()

    def sweep(self, delay_init=1,phase_start=0, phase_stop=2*np.pi, num_phases=50, delay=0.1):
        self.set_phase(0)
        time.sleep(delay_init)
        phase_values = np.linspace(phase_start, phase_stop, num_phases)
        for phase in phase_values:
            print(rf"Phase: {phase/np.pi:.2f}*π", end="\r")
            self.set_phase(float(phase))
            time.sleep(delay)

    def run_optical_sweep(self, powermeters, start_v=0, stop_v=6.1, end_v=0.1, delay=0.02):
        measured_voltage_arr = []
        measured_optical_power_arr = [[] for _ in powermeters] 
        for voltage in np.arange(start_v, stop_v, end_v):
            print('%.2f' % (voltage), end='\r')
            self.chip.qontrol.set_v(self.channel, float(voltage))
            time.sleep(delay)
            measured_voltage_arr.append(self.chip.qontrol.read_v(self.channel))
            for powermeter in range(len(powermeters)):
                measured_optical_power_arr[powermeter].append(powermeters[powermeter].read())     
        self.opt_voltage_arr = np.asarray(measured_voltage_arr)
        self.opt_power_arr = [np.asarray(measured_optical_power_arr[powermeter]) for powermeter in range(len(powermeters))]
        self.visibility = visibility(self.opt_power_arr[0])
        
    def run_optical_fit(self):

        # interpolate the data over the array self.opt_voltage_arr in order to have equal spacing between 1000 points
        # initially, the data is not evely spaced in power, giving more weight to the lower power to fit the curve
        opt_voltage_arr_interp = np.linspace(self.opt_voltage_arr[0], self.opt_voltage_arr[-1], 1000)
        opt_power_arr_interp = np.interp(opt_voltage_arr_interp, self.opt_voltage_arr, self.opt_power_arr[0]) 

        # fit the interpolated data
        popt, _ = fit_optical_curve(
            opt_voltage_arr_interp * (self.a * opt_voltage_arr_interp**2 + self.b * opt_voltage_arr_interp + self.c), 
            opt_power_arr_interp,
            p0=[
                (np.max(opt_power_arr_interp)-np.min(opt_power_arr_interp))/2,
                np.min(opt_power_arr_interp),
                0.02, 
                0]
            )
        
        self.A, self.B, self.omega, self.phi_0 = popt
        self.rho0, self.rho1, self.rho2 = self.omega * self.c, self.omega * self.b, self.omega * self.a

        if self.phase_flip:
            self.phi_0 += np.pi


        # replace data with the interpolated data
        # self.opt_voltage_arr = opt_voltage_arr_interp
        # self.opt_power_arr[0] = opt_power_arr_interp
        
    def find_important_voltages(self):
        self.volt_0    = find_volt(0, self.rho0, self.rho1, self.rho2, self.phi_0)
        self.volt_pi   = find_volt(np.pi, self.rho0, self.rho1, self.rho2, self.phi_0)
        self.volt_pi2  = find_volt(np.pi/2, self.rho0, self.rho1, self.rho2, self.phi_0)
        self.volt_3pi2 = find_volt(3*(np.pi/2), self.rho0, self.rho1, self.rho2, self.phi_0)
    
    def plot_optical_power_voltage(self):
        fig, ax1 = plt.subplots()
        ax1.set_xlabel('V(V)')
        ax1.set_ylabel('Optical Power (μW)')
        for optical_power in self.opt_power_arr:
            ax1.scatter(self.opt_voltage_arr, optical_power*1e6, s=15, c='orchid', label = 'data')
        ax1.tick_params(axis='y')
        ax1.grid(linestyle=':', alpha=0.5)
        ax2 = ax1.twiny()
        ax2.set_xticks([self.volt_0, self.volt_pi, self.volt_pi2, self.volt_3pi2])
        ax2.set_xticklabels(['0', r'$\pi$', r'$\pi/2$', r'$3\pi/2$'])
        ax2.set_xlim(ax1.get_xlim())
        ax2.set_ylabel('Phase (rad)')
        ax2.tick_params(axis='y')
        ax2.grid(color='k', linestyle=':')
        plt.title('MZI '+ self.name)
        plt.savefig(self.folder+f'MZI_{self.name}_optical_voltage.png', dpi=300)
        plt.close()

    def plot_optical_fit(self):
        power = self.opt_voltage_arr * (self.a * self.opt_voltage_arr**2 + self.b * self.opt_voltage_arr + self.c)
        fig, ax1 = plt.subplots()
        ax1.set_xlabel('P(mW)')
        ax1.set_ylabel('Optical Power (μW)')
        ax1.scatter(power, self.opt_power_arr[0]*1e6, s=15, c='orchid', label = 'data')
        ax1.plot(np.linspace(np.min(power), np.max(power), 1000), fit_function(np.linspace(np.min(power), np.max(power), 1000), self.A, self.B, self.omega, self.phi_0)*1e6, 'mediumblue', linewidth=0.7, label='fit')
        ax1.tick_params(axis='y')
        ax1.grid(linestyle=':', alpha=0.5)

        ax2 = ax1.twiny()
        ax2.set_xticks([x * (self.a* x**2 + self.b * x + self.c) for x in [self.volt_0, self.volt_pi, self.volt_pi2, self.volt_3pi2]])
        ax2.set_xticklabels(['0', r'$\pi$', r'$\pi/2$', r'$3\pi/2$'])
        ax2.set_xlim(ax1.get_xlim())
        ax2.set_ylabel('Phase (rad)', color='black')
        ax2.tick_params(axis='y', labelcolor='black')
        ax2.grid(color='k', linestyle=':')

        plt.title('MZI '+ self.name + (" (flipped)" if self.phase_flip else " (normal)"))
        plt.savefig(self.folder+f'MZI_{self.name}_optical.png', dpi=300)
        plt.show()
        plt.close()
    
    def set_phase(self, phi):
        self.chip.qontrol.set_v(self.channel, find_volt(phi, self.rho0, self.rho1, self.rho2, self.phi_0))
        #print(f'Setting phase {phi:.2f} rad on MZI {self.name} (V = {find_volt(phi, self.rho0, self.rho1, self.rho2, self.phi_0):.2f} V)')


 
       
def fit_function(pwr, A, B, omega, phi_0):
    return ((A+B) - A * np.cos(omega * pwr - phi_0))

def fit_optical_curve(pwr, opt_pwr, p0, f=fit_function): #p0=[1e-3, 1e-8, 0.4, 0]
    popt, pcov = curve_fit(f, pwr, opt_pwr, p0=p0)
    if popt[0] < 0:
        p0[3] += np.pi
        popt, pcov = curve_fit(f, pwr, opt_pwr,p0)  
    return popt, pcov
    
def voltage_to_phase(v, phi_0, omega, a, b, c):
    return (omega * v * (a * v**2 + b * v + c) - phi_0)

def compute_roots (rho0, rho1, rho2, P):
    sol1 = -rho1/(3*rho2) - (-3*rho0/rho2 + rho1**2/rho2**2)/(3*(27*P/(2*rho2) - 9*rho0*rho1/(2*rho2**2) + rho1**3/rho2**3 + cmath.sqrt(-4*(-3*rho0/rho2 + rho1**2/rho2**2)**3 + (27*P/rho2 - 9*rho0*rho1/rho2**2 + 2*rho1**3/rho2**3)**2)/2)**(1/3)) - (27*P/(2*rho2) - 9*rho0*rho1/(2*rho2**2) + rho1**3/rho2**3 + cmath.sqrt(-4*(-3*rho0/rho2 + rho1**2/rho2**2)**3 + (27*P/rho2 - 9*rho0*rho1/rho2**2 + 2*rho1**3/rho2**3)**2)/2)**(1/3)/3
    sol2 = -rho1/(3*rho2) - (-3*rho0/rho2 + rho1**2/rho2**2)/(3*(-1/2 - cmath.sqrt(3)*1j/2)*(27*P/(2*rho2) - 9*rho0*rho1/(2*rho2**2) + rho1**3/rho2**3 + cmath.sqrt(-4*(-3*rho0/rho2 + rho1**2/rho2**2)**3 + (27*P/rho2 - 9*rho0*rho1/rho2**2 + 2*rho1**3/rho2**3)**2)/2)**(1/3)) - (-1/2 - cmath.sqrt(3)*1j/2)*(27*P/(2*rho2) - 9*rho0*rho1/(2*rho2**2) + rho1**3/rho2**3 + cmath.sqrt(-4*(-3*rho0/rho2 + rho1**2/rho2**2)**3 + (27*P/rho2 - 9*rho0*rho1/rho2**2 + 2*rho1**3/rho2**3)**2)/2)**(1/3)/3
    sol3 = -rho1/(3*rho2) - (-3*rho0/rho2 + rho1**2/rho2**2)/(3*(-1/2 + cmath.sqrt(3)*1j/2)*(27*P/(2*rho2) - 9*rho0*rho1/(2*rho2**2) + rho1**3/rho2**3 + cmath.sqrt(-4*(-3*rho0/rho2 + rho1**2/rho2**2)**3 + (27*P/rho2 - 9*rho0*rho1/rho2**2 + 2*rho1**3/rho2**3)**2)/2)**(1/3)) - (-1/2 + cmath.sqrt(3)*1j/2)*(27*P/(2*rho2) - 9*rho0*rho1/(2*rho2**2) + rho1**3/rho2**3 + cmath.sqrt(-4*(-3*rho0/rho2 + rho1**2/rho2**2)**3 + (27*P/rho2 - 9*rho0*rho1/rho2**2 + 2*rho1**3/rho2**3)**2)/2)**(1/3)/3
    return [sol1, sol2, sol3]

def find_volt(phi, rho0, rho1, rho2, phi_0):
    Vmax = 24
    phi = phi % (2 * np.pi)

    VoltList = []
    nb_iter = 4 # try to find roots in x intervals below and above the 0-2pi range
    sign = 1 # to scan consecutively the positive and negative ranges
    P = - phi_0 - phi
    
    for iter in range(nb_iter*2+1):
         P +=  iter * sign * (2*np.pi)
         roots = compute_roots(rho0, rho1, rho2, P)
         sign *= -1
         #Choose values that are real and within the accepted v_max range
         for tempVolt in roots:
                 if 0 <= tempVolt.real <= Vmax and abs(tempVolt.imag) < 1e-3:
                     VoltList.append(tempVolt.real)
         # Returns the minimum voltage required
         result = min(VoltList) if VoltList else None
         
    return result
