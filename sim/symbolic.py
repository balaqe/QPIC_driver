import numpy as np
import sympy
import matplotlib.pyplot as plt
from simulator import Clements 

chip = Clements(6)

chip.configure_phase(layer_num=0, shifter=0, phi=0)
chip.configure_phase(layer_num=0, shifter=1, phi=0)
chip.configure_phase(layer_num=0, shifter=2, phi=0)
chip.configure_phase(layer_num=0, shifter=3, phi=0)
chip.configure_phase(layer_num=0, shifter=4, phi=0)
chip.configure_phase(layer_num=0, shifter=5, phi=0)

# chip.configure_mzi(layer_num=1, mzi=0, phi=[np.pi/2, -np.pi/2])
chip.configure_mzi(layer_num=1, mzi=0, phi=[np.pi, -np.pi])
chip.configure_mzi(layer_num=1, mzi=1, phi=[np.pi, -np.pi])
chip.configure_mzi(layer_num=1, mzi=2, phi=[np.pi/2, 3*np.pi/2])

chip.configure_mzi(layer_num=2, mzi=0, phi=[np.pi, -np.pi])
chip.configure_mzi(layer_num=2, mzi=1, phi=[np.pi, np.pi])

# chip.configure_mzi(layer_num=3, mzi=0, phi=[np.pi/2, -np.pi/2])
chip.configure_mzi(layer_num=3, mzi=0, phi=[np.pi, -np.pi])
chip.configure_mzi(layer_num=3, mzi=1, phi=[np.pi, np.pi])
chip.configure_mzi(layer_num=3, mzi=2, phi=[np.pi, -np.pi])

chip.configure_mzi(layer_num=4, mzi=0, phi=[2*np.pi, 0])
chip.configure_mzi(layer_num=4, mzi=1, phi=[np.pi, -np.pi])

# chip.configure_mzi(layer_num=5, mzi=0, phi=[np.pi/2, -np.pi/2])
chip.configure_mzi(layer_num=5, mzi=0, phi=[np.pi, -np.pi])
chip.configure_mzi(layer_num=5, mzi=1, phi=[np.pi, -np.pi])
chip.configure_mzi(layer_num=5, mzi=2, phi=[np.pi, -np.pi])
# chip.configure_mzi(layer_num=5, mzi=2, phi=[np.pi/2, -np.pi/2])

# chip.configure_mzi(layer_num=6, mzi=0, phi=[np.pi/2, -np.pi/2])
chip.configure_mzi(layer_num=6, mzi=0, phi=[0, 0])
chip.configure_mzi(layer_num=6, mzi=1, phi=[np.pi, -np.pi])
# chip.configure_mzi(layer_num=6, mzi=1, phi=[np.pi/2, -np.pi/2])

chip.build()

c = []

for i in range(6):
    c.append(sympy.Symbol(f'c{i}'))

state = chip.execute([c[0], c[1], c[2], c[3], c[4], c[5]])


for i in range(len(state)):
    for a in sympy.preorder_traversal(state[i][0]):
        if isinstance(a, sympy.Float):
            state[i][0] = state[i][0].subs(a, round(a, 1))

print("Resulting state:\n")
print(state)

U = np.array([
    [1, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0],
    [0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 1, 0],
])

print("\n\nUnitary effect:\n")

res = U @ state
print(res)
