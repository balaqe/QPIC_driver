from abc import ABC, abstractmethod
from mpmath.libmp.libintmath import MAX_EULER_CACHE
import numpy as np
# import sympy as sp
from copy import deepcopy

class Mesh_element(ABC):
    @abstractmethod
    def __init__(self, isactive: bool) -> None:
        self.active = False
        self.layer_num: int # Index of the accompanying layer
        self.index: int # Index within a given layer

    @abstractmethod
    def add_phase_offset(self, offset: complex):
        pass

class Mzi(Mesh_element):
    def __init__(self, isactive=True, delta: complex=None, sigma: complex=None):
        self.active = isactive # Switch for indicating whether the mzi is active in the grid

        # self.layer_num: int # Index of the accompanying layer
        # self.index: int # Index within a given layer
        self.upper_pred: Mzi | Phase_shifter # Reference to the upper predecessor
        self.lower_pred: Mzi | Phase_shifter # Lower predecessor
        self.upper_succ: Mzi | Phase_shifter # Upper successor
        self.lower_succ: Mzi | Phase_shifter # Lower successor

        self.sigma = sigma if sigma else 0
        self.delta = delta if delta else 0
        self.phi = [0, 0] # Phi1 and phi2
        self.phase_diff: complex = 0 # Phase difference to be added to lower_pred

        self.get_phase()

    def add_phase_offset(self, offset: complex):
        self.phi[0] += offset
        self.phi[1] += offset

    def get_phase(self):
        self.phi[0] = self.sigma + self.delta
        self.phi[1] = self.sigma - self.delta

class Phase_shifter(Mesh_element):
    def __init__(self, isactive=True):
        self.active = isactive
        #
        # self.layer_num: int # Index of the accompanying layer
        # self.index: int # Index within a given layer
        self.pred: Mzi | Phase_shifter # Reference to the predecessor
        self.succ: Mzi | Phase_shifter # Successor

        self.phi: complex = 0

    def add_phase_offset(self, offset: complex):
        self.phi += offset


class Compiler:
    def __init__(self):
        self.mzi = [] # List of MZIs
        self.phase_shifter = [] # List of phase shifters
        self.og_unitary = None # Original unitary before decomposition

    def compile(self, U):
        self.og_unitary = deepcopy(U)

        x = 0
        y = 0

        x_shape = U.shape[1]
        y_shape = U.shape[0]

        # Initialize MZI structure given the mesh size
        # self.mzi = [[Mzi() for _ in range(x_shape+1)] for _ in range(y_shape)] # x_shape + 1 to account for initial phase shifters
        self.mzi = [[Mzi(isactive=False) for _ in range(x_shape)] for _ in range(y_shape+1)] # y_shape + 1 to account for initial phase shifters
        # self.mesh_elements = [[None for _ in range(x_shape)] for _ in range(y_shape+1)] # y_shape + 1 to account for initial phase shifters
        # self.phase_shifter = [[Phase_shifter() for _ in range(x_shape+1)] for _ in range(y_shape)]
        self.phase_shifter = [[Phase_shifter() for _ in range(x_shape)] for _ in range(y_shape+1)]

        self.mesh_elements = np.empty((y_shape+1, x_shape), dtype=object)

        for elem_num in range(x_shape):
            ps = Phase_shifter()
            ps.layer_num = 0
            ps.index = elem_num
            self.mesh_elements[0][elem_num] = ps
        for layer_num in range(2, y_shape, 2):
            ps1 = Phase_shifter()
            ps1.layer_num = layer_num
            ps1.index = 0
            ps2 = Phase_shifter()
            ps2.layer_num = layer_num
            ps2.index = y_shape-1
            self.mesh_elements[layer_num][0] = ps1
            self.mesh_elements[layer_num][x_shape-1] = ps2
        for layer_num in range(y_shape+1):
            for elem_num in range(x_shape):
                if not self.mesh_elements[layer_num][elem_num]:
                    self.mesh_elements[layer_num][elem_num] = Mzi(isactive=False)

        # V = deepcopy(U)
        V = U.copy().astype(np.complex128)

        for diag_id in range(x_shape-1): # Diagonals starting from the bottom left
            for elem_id in range(diag_id + 1):
            # for elem_id in range(diag_id):
                print(f"diag_id = {diag_id}; elem_id = {elem_id}; element{diag_id + elem_id}({y}, {x})")
                print(f"V = {np.round(V, 2)}")

                # y = y_shape - diag_id + elem_id - 1
                # x = elem_id

                if diag_id % 2 == 0: # j = 1, 3, 5, ... (j starts at 1 while diag_id starts at 0)
                    x = diag_id - elem_id
                    y = y_shape-1 - elem_id
                    print(f"VM")
                    # VM
                    delta = 0
                    if V[y][x] != 0:
                        delta = np.atan(-V[y][x+1] / V[y][x])
                    else:
                        delta = np.pi/2 # Pi/2 shift which produces identity

                    print(f"delta = {delta}")

                    # sigma = np.angle(V[y][x])
                    sigma = 0
                    if y < y_shape-1:
                        sigma = np.angle(V[y+1][x]) - np.angle(V[y][x])
                    else:
                        sigma = -np.angle(V[y][x])

                    phase_diff = np.angle(V[y][x]) - np.angle(V[x][y])

                    # phase_diff = np.angle(V[y][x]) - np.angle(V[y_shape-1 - y][x_shape-1 - x])

                    m = np.array([[np.exp(sigma*1j) * np.sin(delta), np.exp(sigma*1j) * np.cos(delta)],\
                                [np.exp(sigma*1j) * np.cos(delta), -np.exp(sigma*1j) * np.sin(delta)]])

                    print(f"m = {m}")

                    x_offset = x if x < x_shape-2 else x_shape-2
                    y_offset = x_offset

                    M = np.identity(x_shape, dtype=np.complex128)
                    M[y_offset:y_offset+m.shape[0], x_offset:x_offset+m.shape[1]] = m

                    print("Before:")
                    print(np.round(V, 2))
                    V[x][y] *= np.exp(1j*phase_diff)
                    print("After:")
                    print(np.round(V, 2))

                    print(f"M = {np.round(M, 2)}")

                    V = V @ M

                    mzi = Mzi(delta=delta, sigma=sigma)
                    mzi.layer_num = elem_id + 1 # Odd side (start at 1 because layer 0 is phase shifters)
                    mzi.index = diag_id - (mzi.layer_num-1) # Don't divide by 2 to allow half-step offset between layers
                    print(f"MZI added at layer {mzi.layer_num} index {mzi.index}")
                    mzi.phase_diff += phase_diff
                    mzi.active = True
                    self.insert_mzi(mzi)

                else: # j = 0, 2, 4, ...
                    y = y_shape - diag_id + elem_id - 1
                    x = elem_id

                    print(f"MV")
                    # MV
                    delta = 0
                    if V[y][x] != 0:
                        delta = np.atan(V[y-1][x] / V[y][x])
                    else:
                        delta = np.pi/2
                    print(f"delta = {delta}")

                    # sigma = np.angle(V[y][x])
                    sigma = 0
                    if x > 0:
                        sigma = np.angle(V[y][x-1]) - np.angle(V[y][x])
                    else:
                        sigma = -np.angle(V[y][x])

                    phase_diff = np.angle(V[y][x]) - np.angle(V[x][y])

                    m = np.array([[np.exp(sigma*1j) * np.sin(delta), np.exp(sigma*1j) * np.cos(delta)],\
                                [np.exp(sigma*1j) * np.cos(delta), -np.exp(sigma*1j) * np.sin(delta)]])

                    # x_offset = x if x < x_shape-2 else x_shape-2
                    # x_offset = x + 1 if x + 1 < x_shape-2 else x_shape-3
                    x_offset = x+1
                    y_offset = x_offset

                    M = np.identity(x_shape, dtype=np.complex128)
                    M[y_offset:y_offset+m.shape[0], x_offset:x_offset+m.shape[1]] = m

                    V[x][y] *= np.exp(1j*phase_diff)

                    print(f"x_offset = {x_offset}, y_offset = {y_offset}")
                    print(f"M = {np.round(M, 2)}")


                    V = M @ V

                    mzi = Mzi(delta=delta, sigma=sigma)
                    mzi.layer_num = x_shape - elem_id # Even side (don't subtract 1 because the initial phase shifters make the circuit length x_shape+1)
                    # print(f"x_shape = {x_shape}, elem_id = {elem_id}")
                    # mzi.index = y_shape - 2*diag_id - 1 # Don't divide by 2 to allow half-step offset between layers
                    # mzi.index = elem_id*2 + (1-mzi.layer_num%2)
                    diag_offset = y_shape - (diag_id + 2)
                    mzi.index = diag_offset + elem_id

                    mzi.phase_diff += phase_diff
                    mzi.active = True
                    self.insert_mzi(mzi)
        print(f"V = {np.round(V, 2)}")
        self.relax_mesh()

    def insert_mzi(self, mzi: Mzi):
        layer_num = mzi.layer_num
        index = mzi.index
        mzi.get_phase()
        self.mesh_elements[layer_num][index] = mzi

    def relax_mesh(self): # Absorb phase_diffs
        for i, layer in enumerate(self.mesh_elements):
            for j, element in enumerate(layer):
                if not element.active: continue
                if isinstance(element, Phase_shifter): continue # Phase shifters don't have phase_diff
                if element.phase_diff != 0:
                    for k in range(j, len(self.mesh_elements[0])):
                        if self.mesh_elements[i-1][k].active:
                            self.mesh_elements[i-1][k].add_phase_offset(element.phase_diff)
                        if self.mesh_elements[i][k].active:
                            self.mesh_elements[i][k].add_phase_offset(-element.phase_diff)
                    element.phase_diff = 0


# U = np.array( # CNOT
#     [
#         [1, 0, 0, 0],
#         [0, 1, 0, 0],
#         [0, 0, 0, 1],
#         [0, 0, 1, 0]
#     ]
# )
#
def main():
    U = np.array([
        [1, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 1, 0],
    ])

    compiler = Compiler()
    compiler.compile(U)

    i = 0
    for row in compiler.mesh_elements:
        for element in row:
            # print(f"i = {i}: is active? {mzi.active}")
            i += 1
            if not element.active: continue
            if isinstance(element, Mzi): print(f"mzi({element.layer_num}, {element.index}): phi1 = {element.phi[0]}, phi2 = {element.phi[1]}")
            if isinstance(element, Phase_shifter): print(f"phase_shifter({element.layer_num}, {element.index}): phi = {element.phi}")


if __name__=="__main__":
    main()
