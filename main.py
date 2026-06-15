from decomp.tests.arbitrary_q import *
from sim.simulator import Clements
import numpy as np
import sympy as sp
from scipy.stats import unitary_group
from clements import decompose_clements

def is_unitary(m):
    m = np.matrix(m)
    return np.allclose(np.eye(m.shape[0]), m.H * m)

# U = np.array([
#     [1, 0, 0, 0, 0, 0],
#     [0, 1, 0, 0, 0, 0],
#     [0, 0, 1, 0, 0, 0],
#     [0, 0, 0, 1, 0, 0],
#     [0, 0, 0, 0, 0, 1],
#     [0, 0, 0, 0, 1, 0],
# # ])
# U = np.array([
#     [1, 0, 0, 0],
#     [0, 1, 0, 0],
#     [0, 0, 1/np.sqrt(2), 1/np.sqrt(2)],
#     [0, 0, 1/np.sqrt(2), -1/np.sqrt(2)],
# ], dtype=np.complex128)
# U = np.array([
#     [1/np.sqrt(2), 1/np.sqrt(2), 0, 0],
#     [1/np.sqrt(2), -1/np.sqrt(2), 0, 0],
#     [0, 0, 1, 0],
#     [0, 0, 0, 1],
# ])
# U = np.array([
#     [1, 0, 0, 0],
#     [0, 1/np.sqrt(2), 1/np.sqrt(2), 0],
#     [0, 1/np.sqrt(2), -1/np.sqrt(2), 0],
#     [0, 0, 0, 1],
# ])

chip_generator = Clements(4) 
U = chip_generator.generate_random_unitary(seed=0)
print("Initial unitary:")
print(np.round(U, 2))
print(f"Is unitary: {is_unitary(U)}")

# phis, thetas, alphas = decompose_clements(U)

compiler = Compiler()
remultiplied = compiler.compile(U)

chip = Clements(4)

for row in compiler.mesh_elements:
    for element in row:
        if not element.active: continue
        if isinstance(element, Mzi):
            # print(f"mzi({element.layer_num}, {element.index}): phi1 = {element.phi[0]}, phi2 = {element.phi[1]}")
            index = element.index//2 if element.index%2 == 0 else element.index//2+1
            chip.configure_mzi(layer_num=element.layer_num, mzi=index, phi=[element.phi[0], element.phi[1]])
        if isinstance(element, Phase_shifter):
            # print(f"phase_shifter({element.layer_num}, {element.index}): phi = {element.phi}")
            chip.configure_phase(layer_num=element.layer_num, shifter=element.index, phi=element.phi)

chip.build()
# chip.display()

# ex_state = [0, 1, 0, 0]
# # ex_state = [0, 1, 0, 0, 0, 0]
# print("\nResulting state")
# res = chip.execute(ex_state)
# for i in range(len(ex_state)):
#     res[i] = abs(res[i])**2
# print(np.round(res, 2))
#
# print("\nActual state")
# res = U @ np.array(ex_state).reshape(-1, 1)
# for i in range(len(ex_state)):
#     res[i] = abs(res[i])**2
# print(np.round(res, 2))


sym_state = []
for i in range(4):
    sym_state.append(sp.Symbol(f"c{i}"))

# state = chip.execute(sym_state)
state = remultiplied @ np.array(sym_state).reshape(-1, 1)

print("\nResulting state:")
for i in range(len(state)):
    for a in sp.preorder_traversal(state[i][0]):
        if isinstance(a, sp.Float):
            state[i][0] = state[i][0].subs(a, round(a, 1))

# for i in range(len(state)):
#     state[i] = np.abs(state[i])**2

print(state)

state2 = U @ np.array(sym_state).reshape(-1, 1)
# state2 = remultiplied @ np.array(sym_state).reshape(-1, 1)

print("\nActual state")
for i in range(len(state2)):
    for a in sp.preorder_traversal(state2[i][0]):
        if isinstance(a, sp.Float):
            state2[i][0] = state2[i][0].subs(a, round(a, 1))

# for i in range(len(state2)):
#     state2[i] = np.abs(state2[i])**2

print(state2)

# print("\nRemultiplied:")
# print(np.round(remultiplied, 2))

print("\nstate == state2")
print(state == state2)

print(f"IS_UNITARY: {is_unitary(chip.eff_unitary)}")
