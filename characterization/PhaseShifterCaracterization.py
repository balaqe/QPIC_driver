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
# wl_pump = 1546.12
# wl_idler = 1552.52
# wl_signal = 1539.77
# wl_qdot_pump = 1528.6

pm = ThorlabsPowermeter(address='USB0::0x1313::0x8078::P0028334::INSTR', 
                        wavelength=1550)


############# 2. Initialize the chip and the phase shifter ############
#%% Chip

q.reset()
chip = Chip(qontrol=q, folder=r"C:\Users\FTNK-LocalAdm\QPIC_driver\characterization\data\session0\data_master.json")
print(f"Type of chip: {type(chip)}")
print(hasattr(chip, 'mzi_list'))
ps1 = PhaseShifter(chip, name='H11')
ps1.in_port  = 10
ps1.out_port = 10
ps1.channel  = 11
ps1.temperature = 8.54

# test_ps1.chip_config = {
#     'I1': 1
# }

test_ps2 = PhaseShifter(chip, name='H12')

mzi = MZI(chip=chip, name="MZI0", shifter1=ps1, shifter2=test_ps2, folder=r"C:\Users\FTNK-LocalAdm\QPIC_driver\characterization\data\session0")

chip.add_mzi(mzi)


#%%######### 3. Run the electrical characterization ############
# Parameters of the sweep
start_v = 0.0 # Starting voltage
step_v  = 0.1 # Increment voltage
stop_v  = 12 + step_v # Final voltage

# q.set_v(channel=31, voltage=10)

ps1.run_electrical_sweep(start_v, stop_v, step_v, delay=0.02)
ps1.run_electrical_fit()
ps1.plot_electrical_fit()
ps1.plot_electrical_power()


#%%#########
q.reset()
for i in range(43):
    q.set_v(i, 7)


#%%#########

q.reset()
ch2 = 25
q.set_v(ch2, 7)
# q.set_maxv(ch2, 12)
# q.set_maxi(ch2, 50)

print(f"Max voltage: {q.q.vmax[ch2]}")
print(f"Max current: {q.q.imax[ch2]}")

print(f"Device ID: {type(q.q.device_id)}")


#%%#########
#ps1.channel = 7

start_v = 0.0 # Starting voltage
step_v  = 0.1 # Increment voltage
stop_v  = 6 + step_v # Final voltage

# print(f"Max voltage: {q.q.vmax[ch]}")
# print(f"Max current: {q.q.imax[ch]}")

for ch in range(1, 43):
    voltage = []
    current = []
    for i in np.linspace(start_v, stop_v, 20):
        #q.set_v(ps1.channel, float(i))
        q.set_v(ch, float(i))
        voltage.append(q.read_v(ch))
        current.append(q.read_i(ch))
        print(f"voltage: {voltage[-1]}")
    if max(current) > 3:
        plt.plot(voltage, current)
    else:
        plt.plot(voltage, current, label=str(ch))
plt.legend()
plt.show()


# q.set_v(ps1.channel, float(6))
# while(True):
#     current = q.read_i(ps1.channel)
#     print(f"commanded voltage: {6}, actual voltage: {q.read_v(ps1.channel)}, current: {current}")
#     time.sleep(0.5)


#%%######### 4. Run the optical characterization ############
# Parameters of the sweep
# ps1.channel = 7

start_v = 0.0 # Starting voltage
step_v  = 0.05 # Increment voltage
stop_v  = 12 + step_v # Final voltage
# Optical sweep
q.set_v(ps1.channel, 0)
time.sleep(1)

ps1.run_optical_sweep([pm], start_v, stop_v, step_v, delay=0.2)
print(f'Visibility: {ps1.visibility:.3f}')

#%% Optical power fit
ps1.run_optical_fit() 
print(ps1.A, ps1.B, ps1.omega, ps1.phi_0)   
# Find voltages for phases 0, pi, pi/2 and 3pi/2
ps1.find_important_voltages()
# Plots
ps1.plot_optical_power_voltage()
ps1.plot_optical_fit()


########## 5. Save data and close instruments ############
#%% Save data in a json file
chip.save_parameters()
#%% Switch off all the heaters
q.reset()

#%% Close instruments
q.close()
#%% Disable laser
# laser.disable()
