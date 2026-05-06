import numpy as np
import matplotlib.pyplot as plt
from simulator import Clements 

# chip = Clements(2)
#
# # chip.configure_mzi(layer_num=1, mzi=0, phi=[3*np.pi/2, np.pi/2])
# chip.configure_mzi(layer_num=1, mzi=0, phi=[np.pi, 0])
#
# chip.build()
#
# #state = chip.execute([1, 0, 0, 0, 0, 0])
# state = chip.execute([1, 0])
#
# print("Resulting state:\n")
# print(np.round(state, 2))

x = []
y = []

angles = np.linspace(0, 5*np.pi, 1000)
#
for phi in angles:
    chip = Clements(2)
    chip.configure_mzi(layer_num=1, mzi=0, phi=[phi, 0])
    chip.build()
    state = chip.execute([1, 0])
    print(f"angle = {phi}")
    print("Resulting state:\n")
    print(np.round(state, 2))
    x.append(phi)
    y.append(np.absolute(state[0])**2)

plt.plot(x, y)
plt.show()
