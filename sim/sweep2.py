import numpy as np
import matplotlib.pyplot as plt
from simulator import Clements 

chip = Clements(6)

chip.configure_phase(layer_num=0, shifter=0, phi=0)
chip.configure_phase(layer_num=0, shifter=1, phi=0)
chip.configure_phase(layer_num=0, shifter=2, phi=0)
chip.configure_phase(layer_num=0, shifter=3, phi=0)
chip.configure_phase(layer_num=0, shifter=4, phi=0)
chip.configure_phase(layer_num=0, shifter=5, phi=0)

chip.configure_mzi(layer_num=1, mzi=0, phi=[3*np.pi/2, np.pi/2])
chip.configure_mzi(layer_num=1, mzi=1, phi=[0, 0])
chip.configure_mzi(layer_num=1, mzi=2, phi=[0, 0])

chip.configure_mzi(layer_num=2, mzi=0, phi=[0, 0])
chip.configure_mzi(layer_num=2, mzi=1, phi=[0, 0])

chip.configure_mzi(layer_num=3, mzi=0, phi=[3*np.pi/2, np.pi/2])
chip.configure_mzi(layer_num=3, mzi=1, phi=[0, 0])
chip.configure_mzi(layer_num=3, mzi=2, phi=[0, 0])

chip.configure_mzi(layer_num=4, mzi=0, phi=[0, 0])
chip.configure_mzi(layer_num=4, mzi=1, phi=[0, 0])

chip.configure_mzi(layer_num=5, mzi=0, phi=[3*np.pi/2, np.pi/2])
chip.configure_mzi(layer_num=5, mzi=1, phi=[0, 0])
chip.configure_mzi(layer_num=5, mzi=2, phi=[0, 0])

chip.configure_mzi(layer_num=6, mzi=0, phi=[0, 0])
chip.configure_mzi(layer_num=6, mzi=1, phi=[0, 0])

chip.build()

state = chip.execute([1, 0, 0, 0, 0, 0])
# state = chip.execute([1, 0])

print("Resulting state:\n")
print(np.round(state, 2))

x = []
y = []
angles = np.linspace(0, 5*np.pi, 1000)

for phi in angles:
    chip = Clements(6)

    chip.configure_mzi(layer_num=1, mzi=0, phi=[np.pi/2, -np.pi/2])
    chip.configure_mzi(layer_num=3, mzi=0, phi=[np.pi/2, -np.pi/2])
    chip.configure_mzi(layer_num=5, mzi=0, phi=[phi, 0])
    chip.build()

    state = chip.execute([1, 0, 0, 0, 0, 0])
    x.append(phi)
    y.append(np.absolute(state[0])**2)

plt.plot(x, y)
plt.show()
