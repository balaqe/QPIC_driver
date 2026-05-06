"""
Run the electrical and optical caracterization of phase shifters.
Save the data is the associated folder inside the chip folder.
"""
#%%######### Import packages ############ 
from Functions import *

#%% Testing structure
test_chip = Chip(qontrol="nope", folder=r"")
test_chip.load_params("")
# print(f"Type of chip: {type(test_chip)}")
# print(hasattr(test_chip, 'mzi_list'))
test_ps1 = PhaseShifter(test_chip, name='H4')
test_ps1.in_port  = 4
test_ps1.out_port = 3
test_ps1.channel  = 6
test_ps1.phase_flip = True
test_ps1.temperature = 15

test_ps1.chip_config = {
    'I1': 1
}


test_ps2 = PhaseShifter(test_chip, name='H5')
test_ps2.in_port  = 4
test_ps2.out_port = 3
test_ps2.channel  = 6
test_ps2.temperature = 15

test_ps2.chip_config = {
    'I2': 2
}

test_mzi = MZI(chip=test_chip, name="mzi06", shifter1=test_ps1, shifter2=test_ps2, folder=r"")

test_ps3 = PhaseShifter(test_chip, name='H8') # External heater
test_ps3.in_port = 10
test_ps3.out_port = 14
test_ps3.channel = 21
test_ps3.temperature = 15
test_ps3.chip_config = {
    'I2': 3
}

test_chip.add_mzi(test_mzi)
test_chip.add_ext_phase(test_ps3)

test_chip.save_parameters()


############# 1. Initialize the instruments ############
##%% List instruments
#list_devices()
##%% Qontrol, 4=Receiver, 6=Transmitter, 6=Packaged
#q = Qontrol(COM=6)
#q.reset()
##%% Powermeter
#wl_pump = 1546.12
#wl_idler = 1552.52
#wl_signal = 1539.77
#wl_qdot_pump = 1528.6
#
#pm = ThorlabsPowermeter(address='USB0::0x1313::0x8078::P0017921::0::INSTR', 
#                        wavelength=1550)
#
#
############## 2. Initialize the chip and the phase shifter ############
##%% Chip
#chip = Chip(qontrol=q, 
#            folder = r"/Users/balage/Documents/@UNI/@KU/@PHASE3/@PROJ/characterization/data/session0")
##%% Phase shifter
#ps1 = PhaseShifter(chip, name='')
#ps1.in_port  = 4
#ps1.out_port = 3
#ps1.channel  = 6
#
#ps1.chip_config = {
#    'I1': 1
#}
#ps1.chip_config = {key: value * np.pi for key, value in ps1.chip_config.items()}
#chip.set_config(ps1.chip_config)
#
#
##%%######### 3. Run the electrical characterization ############
## Parameters of the sweep
#start_v = 0.0 # Starting voltage
#step_v  = 0.1 # Increment voltage
#stop_v  = 6 + step_v # Final voltage
#
#ps1.run_electrical_sweep(start_v, stop_v, step_v, delay=0.02)
#ps1.run_electrical_fit()
#ps1.plot_electrical_fit()
#ps1.plot_electrical_power()
#
#
##%%######### 4. Run the optical characterization ############
## Parameters of the sweep
#start_v = 0.0 # Starting voltage
#step_v  = 0.05 # Increment voltage
#stop_v  = 6 + step_v # Final voltage
## Optical sweep
#q.set_v(ps1.channel, 0)
#time.sleep(1)
#
#ps1.run_optical_sweep([pm], start_v, stop_v, step_v, delay=0.2)
#print(f'Visibility: {ps1.visibility:.3f}')
##%% Optical power fit
#ps1.run_optical_fit() 
#print(ps1.A, ps1.B, ps1.omega, ps1.phi_0)   
## Find voltages for phases 0, pi, pi/2 and 3pi/2
#ps1.find_important_voltages()
## Plots
#ps1.plot_optical_power_voltage()
#ps1.plot_optical_fit()
#
#
########### 5. Save data and close instruments ############
##%% Save data in a json file
#ps1.save_parameters()
##%% Switch off all the heaters
#q.reset()
#
##%% Close instruments
#q.close()
##%% Disable laser
## laser.disable()
