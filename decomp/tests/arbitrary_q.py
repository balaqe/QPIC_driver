import numpy as np
import sympy as sp
import deepcopy

class Mzi:
    def __init__(self):
        self.layer_num = None # Index of the accompanying layer
        self.index = None # Index within a given layer
        self.upper_pred = None # Reference to the upper predecessor
        self.lower_pred = None # Lower predecessor
        self.upper_succ = None # Upper successor
        self.lower_succ = None # Lower successor

        self.sigma = None
        self.delta = None
        self.phi = [complex, complex] # Phi1 and phi2
        self.phase_diff = 0 # Phase difference to be added to lower_pred

class Phase_shifter:
    def __init__(self):
        self.layer_num = None # Index of the accompanying layer
        self.index = None # Index within a given layer
        self.pred = None # Reference to the predecessor
        self.succ = None # Successor

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
        self.mzi = [[None for _ in range(x_shape+1)] for _ in range(y_shape)] # x_shape + 1 to account for initial phase shifters
        self.phase_shifter = [[None for _ in range(x_shape+1)] for _ in range(y_shape)]

        V = deepcopy(U)

        for diag_id in range(x_shape): # Diagonals starting from the bottom left
            for elem_id in range(diag_id + 1):
                y = y_shape - diag_id + elem_id - 1
                x = elem_id
                print(f"diag_id = {diag_id}; elem_id = {elem_id}; element{diag_id + elem_id}({x}, {y})")

                if diag_id % 2 == 0: # j = 1, 3, 5, ... (j starts at 1 while diag_id starts at 0)
                    # VM
                    delta = np.atan(V[x][y] / V[x+1][y])
                    sigma = sympy.angle(V[x][y])
                    phase_diff = V[x][y] - V[x_shape-1 - x][y_shape-1 - y]

                    m = np.array([[np.exp(sigma*1j) * np.sin(delta), np.exp(sigma*1j) * np.cos(delta)],\
                                [np.exp(sigma*1j) * np.cos(delta), -np.exp(sigma*1j) * np.sin(delta)]])

                    x_offset = x
                    y_offset = x_offset

                    M = np.identity(x_shape)
                    M[x_offset:x_offset+m.shape[0], y_offset:y_offset+m.shape[1]] = m

                    V = V @ M

                    mzi = Mzi()
                    mzi.layer_num = elem_id + 1 # Odd side (start at 1 because layer 0 is phase shifters)
                    mzi.index = diag_id # Don't divide by 2 to allow half-step offset between layers
                    mzi.delta = delta
                    mzi.sigma = sigma
                    phase_diff = phase_diff
                    self.insert_mzi(mzi)

                else: # j = 0, 2, 4, 
                    # MV

                    mzi = Mzi()
                    mzi.layer_num = x_shape - elem_id - 1 # Even side
                    mzi.index = y_shape - diag_id - 1 # Don't divide by 2 to allow half-step offset between layers
                    mzi.delta = delta
                    mzi.sigma = sigma
                    phase_diff = phase_diff
                    self.insert_mzi(mzi)



    def insert_mzi(self, mzi: Mzi):
        layer_num = mzi.layer_num
        index = mzi.index

        if layer_num > 2: # Connect predecessors
            if index > 0:
                if self.mzi[layer_num-1][index-1] != None:
                    mzi.upper_pred = self.mzi[layer_num-1][index-1]
                    self.mzi[layer_num-1][index-1].lower_succ = mzi
                elif self.phase_shifter[layer_num-1][index-1] != None:
                    mzi.upper_pred = self.phase_shifter[layer_num-1][index-1]
                    self.phase_shifter[layer_num-1][index-1].succ = mzi

            if index < self.mzi.shape[0] - 1:
                if self.mzi[layer_num-1][index+1] != None:
                    mzi.lower_pred = self.mzi[layer_num-1][index+1]
                    self.mzi[layer_num-1][index+1].upper_succ = mzi
                elif self.phase_shifter[layer_num-1][index+1] != None:
                    mzi.upper_pred = self.phase_shifter[layer_num-1][index+1]
                    self.phase_shifter[layer_num-1][index+1].succ = mzi

        if layer_num < self.mzi.shape[1] - 1: # Connect successors
            if index > 0:
                if self.mzi[layer_num+1][index-1] != None:
                    mzi.upper_pred = self.mzi[layer_num+1][index-1]
                    self.mzi[layer_num+1][index-1].lower_succ = mzi
                elif self.phase_shifter[layer_num+1][index-1] != None:
                    mzi.upper_pred = self.phase_shifter[layer_num+1][index-1]
                    self.phase_shifter[layer_num+1][index-1].succ = mzi

            if index < self.mzi.shape[0] - 1:
                if self.mzi[layer_num+1][index+1] != None:
                    mzi.lower_pred = self.mzi[layer_num+1][index+1]
                    self.mzi[layer_num+1][index+1].upper_succ = mzi
                elif self.phase_shifter[layer_num+1][index+1] != None:
                    mzi.upper_pred = self.phase_shifter[layer_num+1][index+1]
                    self.phase_shifter[layer_num+1][index+1].succ = mzi




U = np.array( # CNOT
    [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 1, 0]
    ]
)

compiler = Compiler()
compiler.compile(U)



# x = 0
# y = 0
#
# x_shape = U.shape[1]
# y_shape = U.shape[0]
#
# V = deepcopy(U)
#
# for diag_id in range(x_shape): # Diagonals starting from the bottom left
#     for elem_id in range(diag_id + 1):
#         y = y_shape - diag_id + elem_id - 1
#         x = elem_id
#         print(f"diag_id = {diag_id}; elem_id = {elem_id}; element{diag_id + elem_id}({x}, {y})")
#
#         if diag_id % 2 == 0: # j = 1, 3, 5, ... (j starts at 1 while diag_id starts at 0)
#             # VM
#             delta = np.atan(V[x][y] / V[x+1][y])
#             sigma = sympy.angle(V[x][y])
#             phase_diff = V[x][y] - V[x_shape-1 - x][y_shape-1 - y]
#
#             m = np.array([[np.exp(sigma*1j) * np.sin(delta), np.exp(sigma*1j) * np.cos(delta)],\
#                           [np.exp(sigma*1j) * np.cos(delta), -np.exp(sigma*1j) * np.sin(delta)]])
#
#             x_offset = x
#             y_offset = x_offset
#
#             M = np.identity(x_shape)
#             M[x_offset:x_offset+m.shape[0], y_offset:y_offset+m.shape[1]] = m
#
#             V = V @ M
#
#         else: # j = 0, 2, 4, 
#             # MV
