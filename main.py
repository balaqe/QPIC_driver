from decomp.tests.arbitrary_q import *
from sim.simulator import Clements
import numpy as np
import sympy as sp

# U = np.array([
#     [1, 0, 0, 0, 0, 0],
#     [0, 1, 0, 0, 0, 0],
#     [0, 0, 1, 0, 0, 0],
#     [0, 0, 0, 1, 0, 0],
#     [0, 0, 0, 0, 0, 1],
#     [0, 0, 0, 0, 1, 0],
# ])
U = np.array([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 1, 0],
])

compiler = Compiler()
compiler.compile(U)

chip = Clements(4)

for row in compiler.mesh_elements:
    for element in row:
        if not element.active: continue
        if isinstance(element, Mzi):
            print(f"mzi({element.layer_num}, {element.index}): phi1 = {element.phi[0]}, phi2 = {element.phi[1]}")
            index = element.index//2 if element.index%2 == 0 else element.index//2+1
            chip.configure_mzi(layer_num=element.layer_num, mzi=index, phi=[element.phi[0], element.phi[1]])
        if isinstance(element, Phase_shifter):
            print(f"phase_shifter({element.layer_num}, {element.index}): phi = {element.phi}")
            chip.configure_phase(layer_num=element.layer_num, shifter=element.index, phi=element.phi)

chip.build()
# chip.display()


sym_state = []
for i in range(4):
    sym_state.append(sp.Symbol(f"c{i}"))

state = chip.execute(sym_state)

print("Resulting state:\n")
for i in range(len(state)):
    for a in sp.preorder_traversal(state[i][0]):
        if isinstance(a, sp.Float):
            state[i][0] = state[i][0].subs(a, round(a, 1))
print(state)

state2 = U @ np.array(sym_state).reshape(-1, 1)

print("Actual state")
for a in sp.preorder_traversal(state2[0]):
    if isinstance(a, sp.Float):
        state2[0] = state2[0].subs(a, round(a, 1))
print(state2)

