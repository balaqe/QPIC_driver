"""
Run the electrical and optical caracterization of phase shifters.
Save the data is the associated folder inside the chip folder.
"""
#%%######### Import packages ############ 
from Functions import *


############ 1. Initialize the instruments ############
#%% List instruments
list_devices()
#%% Qontrol, 4=Receiver, 6=Transmitter, 6=Packaged
q = Qontrol(COM=4)
q.reset()
#%% Powermeter
wl_pump = 1546.12
wl_idler = 1552.52
wl_signal = 1539.77
wl_qdot_pump = 1528.6

pm = ThorlabsPowermeter(address='USB0::0x1313::0x8078::P0028334::0::INSTR', 
                       wavelength=1550)


############# 2. Initialize the chip and the phase shifter ############
#%% Chip
chip = Chip(qontrol=q, 
           folder = r"C:\Users\FTNK-LocalAdm\QPIC_driver\characterization\data\session0")

chip.set_config(
    {
        'MZI0': {
            'shifter1': 1,
        },
    }
)

ps = PhaseShifter(chip, name='H16')
ps.in_port  = 10
ps.out_port = 10
ps.channel  = 16
ps.temperature = 15
ps.phase_flip = True

mzi = None
mzi_name = "MZI3"

top_shifter = True
if top_shifter:
    ps2 = PhaseShifter(chip, name='')
    mzi = MZI(chip, name=mzi_name, shifter1=ps, shifter2=ps2)
else:
    mzi = chip.mzi_dict[mzi_name]
    mzi.shifter[1] = ps

chip.add_mzi(mzi)

# ps1.chip_config = {
#    'H1': 1
# }
# ps1.chip_config = {key: value * np.pi for key, value in ps1.chip_config.items()}

# chip.set_config(ps1.chip_config)


#%%######### 3. Run the electrical characterization ############
# Parameters of the sweep
start_v = 0.0 # Starting voltage
step_v  = 0.1 # Increment voltage
stop_v  = 12 + step_v # Final voltage

ps.run_electrical_sweep(start_v, stop_v, step_v, delay=0.02)
ps.run_electrical_fit()
ps.plot_electrical_fit()
ps.plot_electrical_power()


#%%######### 4. Run the optical characterization ############
# Parameters of the sweep
start_v = 0.0 # Starting voltage
step_v  = 0.05 # Increment voltage
stop_v  = 12 + step_v # Final voltage
# Optical sweep
q.set_v(ps.channel, 0)
time.sleep(1)

ps.run_optical_sweep([pm], start_v, stop_v, step_v, delay=0.05)
print(f'Visibility: {ps.visibility:.3f}')


#%% Optical power fit

ps.run_optical_fit() 
print(ps.A, ps.B, ps.omega, ps.phi_0)   
# Find voltages for phases 0, pi, pi/2 and 3pi/2
ps.find_important_voltages()
# Plots
ps.plot_optical_power_voltage()
ps.plot_optical_fit()


########## 5. Save data and close instruments ############
#%% Save data in a json file
# ps1.save_parameters()
chip.save_parameters()
#%% Switch off all the heaters
q.reset()

#%% Close instruments
q.close()
#%% Disable laser
# laser.disable()
