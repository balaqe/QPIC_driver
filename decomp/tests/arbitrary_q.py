import numpy as np
# import sympy as sp
from copy import deepcopy

class Mzi:
    def __init__(self, isactive = True):
        self.active = isactive # Switch for indicating whether the mzi is active in the grid

        self.layer_num: int # Index of the accompanying layer
        self.index: int # Index within a given layer
        self.upper_pred: Mzi | Phase_shifter # Reference to the upper predecessor
        self.lower_pred: Mzi | Phase_shifter # Lower predecessor
        self.upper_succ: Mzi | Phase_shifter # Upper successor
        self.lower_succ: Mzi | Phase_shifter # Lower successor

        self.sigma: complex
        self.delta: complex
        self.phi = [complex, complex] # Phi1 and phi2
        self.phase_diff: complex # Phase difference to be added to lower_pred

    def get_phases(self):
        self.phi[0] = self.sigma + self.delta
        self.phi[1] = self.sigma - self.delta

class Phase_shifter:
    def __init__(self):
        self.active = False

        self.layer_num: int # Index of the accompanying layer
        self.index: int # Index within a given layer
        self.pred: Mzi | Phase_shifter # Reference to the predecessor
        self.succ: Mzi | Phase_shifter # Successor

        self.phi = None


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
        # self.phase_shifter = [[Phase_shifter() for _ in range(x_shape+1)] for _ in range(y_shape)]
        self.phase_shifter = [[Phase_shifter() for _ in range(x_shape)] for _ in range(y_shape+1)]

        V = deepcopy(U)

        for diag_id in range(x_shape-1): # Diagonals starting from the bottom left
            for elem_id in range(diag_id + 1):
            # for elem_id in range(diag_id):
                y = y_shape - diag_id + elem_id - 1
                x = elem_id
                # print(f"diag_id = {diag_id}; elem_id = {elem_id}; element{diag_id + elem_id}({x}, {y})")

                if diag_id % 2 == 0: # j = 1, 3, 5, ... (j starts at 1 while diag_id starts at 0)
                    # VM
                    delta = 0
                    if V[y][x+1] != 0:
                        delta = np.atan(-V[y][x] / V[y][x+1])
                    else:
                        delta = np.pi # Pi shift which produces identity

                    sigma = np.angle(V[y][x])
                    phase_diff = np.angle(V[y][x]) - np.angle(V[x][y])
                    # phase_diff = np.angle(V[y][x]) - np.angle(V[y_shape-1 - y][x_shape-1 - x])

                    m = np.array([[np.exp(sigma*1j) * np.sin(delta), np.exp(sigma*1j) * np.cos(delta)],\
                                [np.exp(sigma*1j) * np.cos(delta), -np.exp(sigma*1j) * np.sin(delta)]])

                    x_offset = x if x < x_shape-2 else x_shape-2
                    y_offset = x_offset

                    M = np.identity(x_shape, dtype="complex")
                    M[y_offset:y_offset+m.shape[0], x_offset:x_offset+m.shape[1]] = m

                    V = V @ M

                    mzi = Mzi()
                    mzi.layer_num = elem_id + 1 # Odd side (start at 1 because layer 0 is phase shifters)
                    mzi.index = diag_id - (mzi.layer_num-1) # Don't divide by 2 to allow half-step offset between layers
                    mzi.delta = delta
                    mzi.sigma = sigma
                    mzi.phase_diff = phase_diff
                    mzi.active = True
                    mzi.get_phases()
                    self.insert_mzi(mzi)

                else: # j = 0, 2, 4, ...
                    # MV
                    delta = 0
                    if V[y][x] != 0:
                        delta = np.atan(-V[y-1][x] / V[y][x])
                    else:
                        delta = np.pi
                    sigma = np.angle(V[y][x])
                    phase_diff = np.angle(V[y][x]) - np.angle(V[x][y])

                    m = np.array([[np.exp(sigma*1j) * np.sin(delta), np.exp(sigma*1j) * np.cos(delta)],\
                                [np.exp(sigma*1j) * np.cos(delta), -np.exp(sigma*1j) * np.sin(delta)]])

                    x_offset = x if x < x_shape-2 else x_shape-2
                    y_offset = x_offset

                    M = np.identity(x_shape, dtype="complex")
                    M[y_offset:y_offset+m.shape[0], x_offset:x_offset+m.shape[1]] = m

                    V = M @ V

                    mzi = Mzi()
                    mzi.layer_num = x_shape - elem_id # Even side (don't subtract 1 because the initial phase shifters make the circuit length x_shape+1)
                    # print(f"x_shape = {x_shape}, elem_id = {elem_id}")
                    # mzi.index = y_shape - 2*diag_id - 1 # Don't divide by 2 to allow half-step offset between layers
                    # mzi.index = elem_id*2 + (1-mzi.layer_num%2)
                    diag_offset = y_shape - (diag_id + 2)
                    # print(f"diag_offset = {diag_offset}")
                    mzi.index = diag_offset + elem_id

                    # print(f"y_shape = {y_shape}, diag_id = {elem_id}")
                    print(f"diag_id = {diag_id}, elem_id = {elem_id}, layer_num = {mzi.layer_num}, index = {mzi.index}")
                    mzi.delta = delta
                    mzi.sigma = sigma
                    mzi.phase_diff = phase_diff
                    mzi.active = True
                    mzi.get_phases()
                    self.insert_mzi(mzi)

        self.relax_phases()


    def insert_mzi(self, mzi: Mzi):
        layer_num = mzi.layer_num
        index = mzi.index
        self.mzi[layer_num][index] = mzi
        # print(f"is active? {self.mzi[layer_num][index].active}")

        if layer_num > 0: # Connect predecessors
            if index > 0:
                if self.mzi[layer_num-1][index-1].active:
                    mzi.upper_pred = self.mzi[layer_num-1][index-1]
                    self.mzi[layer_num-1][index-1].lower_succ = mzi
                elif self.phase_shifter[layer_num-1][index-1].active:
                    mzi.upper_pred = self.phase_shifter[layer_num-1][index-1]
                    self.phase_shifter[layer_num-1][index-1].succ = mzi

            if index < len(self.mzi) - 1: # index < self.mzi.shape[0] - 1
                if self.mzi[layer_num-1][index+1].active:
                    mzi.lower_pred = self.mzi[layer_num-1][index+1]
                    self.mzi[layer_num-1][index+1].upper_succ = mzi
                elif self.phase_shifter[layer_num-1][index+1].active:
                    mzi.upper_pred = self.phase_shifter[layer_num-1][index+1]
                    self.phase_shifter[layer_num-1][index+1].succ = mzi

        if layer_num < len(self.mzi[0]) - 1: # Connect successors (all rows are the same length so we check the first)
            if index > 0:
                # print(f"DEBUG: layer_num = {layer_num}")
                # print(f"DEBUG: len(self.mzi[0] - 1 = {len(self.mzi[0]) -  1}")
                if self.mzi[layer_num+1][index-1].active:
                    mzi.upper_pred = self.mzi[layer_num+1][index-1]
                    self.mzi[layer_num+1][index-1].lower_succ = mzi
                elif self.phase_shifter[layer_num+1][index-1].active:
                    mzi.upper_pred = self.phase_shifter[layer_num+1][index-1]
                    self.phase_shifter[layer_num+1][index-1].succ = mzi

            if index < len(self.mzi) - 1:
                if self.mzi[layer_num+1][index+1].active:
                    mzi.lower_pred = self.mzi[layer_num+1][index+1]
                    self.mzi[layer_num+1][index+1].upper_succ = mzi
                elif self.phase_shifter[layer_num+1][index+1].active:
                    mzi.upper_pred = self.phase_shifter[layer_num+1][index+1]
                    self.phase_shifter[layer_num+1][index+1].succ = mzi

    def relax_phases(self):
        for layer_num in reversed(range(len(self.mzi))):
            for mzi in self.mzi[layer_num]: # select layer
                if not mzi.active or mzi.phase_diff == 0: continue
                phase_diff = mzi.phase_diff
                for index in range(mzi.index+2, len(self.mzi[0]), 2):
                    print(f"index = {index}, len(self.mzi) = {len(self.mzi[0])}")
                    print(f"type of phi[0]: {self.mzi[layer_num][index].phi[0]}, type of phase_diff: {}")
                    self.mzi[layer_num][index].phi[0] -= mzi.phase_diff
                    self.mzi[layer_num][index].phi[1] -= mzi.phase_diff
                for index in range(0, len(self.mzi), 2):
                    self.mzi[layer_num-1][index].phi[0] += mzi.phase_diff
                    self.mzi[layer_num-1][index].phi[1] += mzi.phase_diff




# U = np.array( # CNOT
#     [
#         [1, 0, 0, 0],
#         [0, 1, 0, 0],
#         [0, 0, 0, 1],
#         [0, 0, 1, 0]
#     ]
# )
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
for row in compiler.mzi:
    for mzi in row:
        # print(f"i = {i}: is active? {mzi.active}")
        i += 1
        if not mzi.active: continue
        print(f"mzi({mzi.layer_num}, {mzi.index}): phi1 = {mzi.phi[0]}, phi2 = {mzi.phi[1]}")
        if mzi.phase_diff != 0:
            print(f"WARNING! phase diff is not zero on mzi({mzi.layer_num}, {mzi.index}).phase_diff = {mzi.phase_diff}")
