from abc import ABC, abstractmethod
from mpmath.libmp.libintmath import MAX_EULER_CACHE
import numpy as np
from numpy import complex128, isin, linalg
# import sympy as sp
from copy import deepcopy

class Mesh_element(ABC):
    @abstractmethod
    def __init__(self) -> None:
        self.layer_num: int # Index of the accompanying layer
        self.index: int # Index within a given layer

    @abstractmethod
    def add_phase_offset(self, offset: np.complex128):
        pass

class Mzi(Mesh_element):
    def __init__(self, delta: np.complex128=None, sigma: np.complex128=None):
        # self.layer_num: int # Index of the accompanying layer
        # self.index: int # Index within a given layer

        self.sigma = sigma if sigma else 0
        self.delta = delta if delta else 0
        self.phi = [0, 0] # Phi1 and phi2
        self.phase_diff: complex = 0 # Phase difference to be added to lower_pred

        self.get_phase()

    def add_phase_offset(self, offset: np.complex128):
        self.phi[0] += offset
        self.phi[1] += offset

    def get_phase(self):
        self.phi[0] = self.sigma + self.delta
        self.phi[1] = self.sigma - self.delta

class Phase_shifter(Mesh_element):
    def __init__(self):
        # self.layer_num: int # Index of the accompanying layer
        # self.index: int # Index within a given layer

        self.phi: complex = 0

    def add_phase_offset(self, offset: complex):
        self.phi += offset


class Compiler:
    def __init__(self):
        self.mzi = [] # List of MZIs
        self.phase_shifter = [] # List of phase shifters
        self.og_unitary = None # Original unitary before decomposition
        self.phantom_phases = []

    def compile(self, U):
        self.og_unitary = deepcopy(U)

        x = 0
        y = 0

        x_shape = U.shape[1]
        y_shape = U.shape[0]
        self.num_layers = x_shape*2

        T = []
        T_dagger = []
        bruh = np.identity(x_shape, dtype=np.complex128)
        bruh_dagger = np.identity(x_shape, dtype=np.complex128)

        self.mesh_elements = np.empty((self.num_layers, y_shape), dtype=object)

        V = np.conjugate(U.copy().astype(np.complex128))
        # V = -1*(U.copy().astype(np.complex128))

        for diag_id in range(x_shape-1): # Diagonals starting from the bottom left
            for elem_id in range(diag_id + 1):
                if diag_id % 2 == 0: # j = 1, 3, 5, ... (j starts at 1 while diag_id starts at 0)
                    # VM
                    x = diag_id - elem_id
                    y = y_shape-1 - elem_id

                    sigma = 0
                    # if y < y_shape-1:
                    #     sigma = np.angle(V[y][x+1]) - np.angle(V[y][x])
                    # else:
                    #     sigma = -np.angle(V[y][x])

                    phase_diff = 0
                    if x < x_shape-1: phase_diff = np.angle(V[y][x+1]) - np.angle(V[y][x])
                    else: phase_diff = -np.angle(V[y][x])

                    P = np.identity(x_shape, dtype=np.complex128)
                    P[x][x] = np.exp(1j*phase_diff)
                    # print("P:")
                    # print(np.round(P, 2))
                    V = V @ P
                    T_dagger.append(P)
                    bruh_dagger = bruh_dagger @ P

                    delta = 0
                    if V[y][x] != 0:
                        delta = np.atan(-V[y][x+1] / V[y][x], dtype=np.complex128)
                    else:
                        delta = np.pi/2 # Pi/2 shift which produces identity

                    m = np.array([[np.exp(sigma*1j) * np.sin(delta), np.exp(sigma*1j) * np.cos(delta)],\
                                [np.exp(sigma*1j) * np.cos(delta), -np.exp(sigma*1j) * np.sin(delta)]], dtype=np.complex128)

                    x_offset = x if x < x_shape-2 else x_shape-2
                    y_offset = x_offset

                    M = np.identity(x_shape, dtype=np.complex128)
                    M[y_offset:y_offset+m.shape[0], x_offset:x_offset+m.shape[1]] = m

                    # print("M:")
                    # print(np.round(M, 2))

                    V = V @ M
                    T_dagger.append(M)
                    bruh_dagger = bruh_dagger @ M

                    mzi = Mzi(delta=delta, sigma=sigma)
                    mzi.layer_num = elem_id*2 + 1 # Odd side (start at 1 because layer 0 is phase shifters)
                    mzi.index = diag_id - (mzi.layer_num-1) // 2 # Don't divide by 2 to allow half-step offset between layers
                    mzi.phase_diff += phase_diff
                    self.insert_mzi(mzi)

                    print(f"MZI added at layer {mzi.layer_num} and index {mzi.index}")

                    phantom_shifter = Phase_shifter()
                    phantom_shifter.phi = phase_diff
                    phantom_shifter.index = mzi.index
                    phantom_shifter.layer_num = mzi.layer_num-1
                    self.phantom_phases.append(phantom_shifter)
                    self.mesh_elements[phantom_shifter.layer_num][phantom_shifter.index] = phantom_shifter
                    print(f"phantom phase added at layer {phantom_shifter.layer_num} and index {phantom_shifter.index}")

                else: # j = 2, 4, ...
                    # MV
                    y = y_shape - diag_id + elem_id - 1
                    x = elem_id

                    sigma = 0
                    # if x > 0:
                    #     sigma = np.angle(V[y][x-1]) - np.angle(V[y][x])
                    # else:
                    #     sigma = -np.angle(V[y][x])
                    #
                    phase_diff = 0
                    if y > 0: phase_diff = np.angle(V[y-1][x]) - np.angle(V[y][x])
                    else: phase_diff = -np.angle(V[y][x])

                    P = np.identity(x_shape, dtype=np.complex128)
                    P[y][y] = np.exp(1j*phase_diff)
                    # print("P:")
                    # print(np.round(P, 2))
                    V = P @ V
                    T.append(P)
                    bruh = P @ bruh

                    delta = 0
                    if V[y][x] != 0:
                        delta = np.atan(V[y-1][x] / V[y][x], dtype=np.complex128)
                    else:
                        delta = np.pi/2

                    m = np.array([[np.exp(sigma*1j) * np.sin(delta), np.exp(sigma*1j) * np.cos(delta)],\
                                [np.exp(sigma*1j) * np.cos(delta), -np.exp(sigma*1j) * np.sin(delta)]], dtype=np.complex128)

                    x_offset = x_shape-diag_id-2 +x
                    y_offset = x_offset

                    M = np.identity(x_shape, dtype=np.complex128)
                    M[y_offset:y_offset+m.shape[0], x_offset:x_offset+m.shape[1]] = m
                    # print("M:")
                    # print(np.round(M, 2))

                    V = M @ V
                    T.append(M)
                    bruh = M @ bruh

                    mzi = Mzi(delta=delta, sigma=sigma)
                    mzi.layer_num = self.num_layers-1 - elem_id*2 # Even side (don't subtract 1 because the initial phase shifters make the circuit length x_shape+1)
                    diag_offset = y_shape - (diag_id + 2)
                    mzi.index = diag_offset + elem_id

                    print(f"MZI added at layer {mzi.layer_num} and index {mzi.index}")

                    mzi.phase_diff += phase_diff
                    self.insert_mzi(mzi)

                    # print(f"MZI added at layer {mzi.layer_num} and index {mzi.index}")

                    phantom_shifter = Phase_shifter()
                    phantom_shifter.phi = phase_diff
                    phantom_shifter.index = mzi.index
                    phantom_shifter.layer_num = mzi.layer_num-1
                    self.phantom_phases.append(phantom_shifter)

                    self.mesh_elements[phantom_shifter.layer_num][phantom_shifter.index] = phantom_shifter

                    print(f"phantom phase added at layer {phantom_shifter.layer_num} and index {phantom_shifter.index}")

        # Get alphas
        # if self.og_unitary.shape[0] % 2 == 0: # Middle diagonal of mesh goes back down
        #     target_layer = self.num_layers-1
        #     for i in range(V.shape[0]):
        #         if isinstance(self.mesh_elements[target_layer][i], Mzi):
        #             target_layer -= 1
        #         phi = np.angle(V[i][i])
        #         if self.mesh_elements[target_layer][i]:
        #             self.mesh_elements[target_layer][i].add_phase_offset(phi)
        #             print(f"Alpha phase shifter added on layer {target_layer} and index {i}")
        #         else:
        #             phantom_shifter = Phase_shifter()
        #             phantom_shifter.phi = phi
        #             phantom_shifter.index = i
        #             phantom_shifter.layer_num = target_layer
        #             self.phantom_phases.append(phantom_shifter)
        #             self.mesh_elements[phantom_shifter.layer_num][phantom_shifter.index] = phantom_shifter
        #             print(f"Alpha phase shifter added on layer {phantom_shifter.layer_num} and index {phantom_shifter.index}")
        #         target_layer -= 1
        # else: # Middle diagonal of mesh goes forward up
        #     target_layer = 1
        #     for i in range(V.shape[0]):
        #         if isinstance(self.mesh_elements[target_layer][self.og_unitary.shape[0]-1 - i], Mzi):
        #             target_layer += 1
        #         phantom_shifter = Phase_shifter()
        #         phantom_shifter.phi = np.angle(V[i][i])
        #         phantom_shifter.index = self.og_unitary.shape[0]-1 - i
        #         phantom_shifter.layer_num = target_layer
        #         self.phantom_phases.append(phantom_shifter)
        #         self.mesh_elements[phantom_shifter.layer_num][phantom_shifter.index] = phantom_shifter
        #         print(f"Alpha phase shifter added on layer {phantom_shifter.layer_num} and index {phantom_shifter.index}")
        #         target_layer += 1

        #     # self.mesh_elements[0][i].add_phase_offset(phase)
        #
        print("V:")
        print(f"{np.round(V, 2)}")
        # self.relax_mesh()

        pre = np.identity(x_shape, dtype=np.complex128)
        post = np.identity(x_shape, dtype=np.complex128)
        # print("Saved elements")

        # i = 0
        # while i < len(T):
        #     p = T[i]
        #     m = T[i+1]
        #     print(f"p:")
        #     print(np.round(p, 2))
        #     print(f"m:")
        #     print(np.round(m, 2))
        #     print("corrected:")
        #     print(np.round(p @ m, 2))
        #     i += 2
        #
        print("REMULTIPLY ORDER")
        for t in T:
            print("T[i]")
            print(np.round(t, 2))
        for t_dagger in reversed(T_dagger):
            print("T_dagger[i]")
            print(np.round(t_dagger, 2))

        for t in T:
            pre =  pre @ t
            # print("T[i]")
            # print(np.round(t, 2))


        # j = 0
        # while j < len(T_dagger):
        #     p = T_dagger[j]
        #     m = T_dagger[j+1]
        #     print(f"p:")
        #     print(np.round(p, 2))
        #     print(f"m:")
        #     print(np.round(m, 2))
        #     print("corrected:")
        #     print(np.round(p @ m, 2))
        #     j += 2

        for t_dagger in T_dagger:
            post = t_dagger @ post
            # print("T_dagger[i]")
            # print(np.round(t_dagger, 2))

        # inv_pre = linalg.inv(bruh)
        # inv_post = linalg.inv(bruh_dagger)

        # print("\ninv_pre:")
        # print(np.round(inv_pre, 2))
        # print("\ninv_post:")
        # print(np.round(inv_post, 2))
        # print("\npre:")
        # print(np.round(pre, 2))
        # print("\npost:")
        # print(np.round(post, 2))

        # mid_P = np.identity(x_shape, dtype=np.complex128)
        # for i in range(x_shape):
        #     angle = np.angle(V[i][i])
        #     mid_P[i][i] = np.exp(-1j*angle)

        # return inv_pre @ mid_P @ inv_post
        # return pre @ mid_P @ post
        return pre @ post

    def insert_mzi(self, mzi: Mzi):
        layer_num = mzi.layer_num
        index = mzi.index
        mzi.get_phase()
        self.mesh_elements[layer_num][index] = mzi

    def print_elements(self):
        for layer in self.mesh_elements:
            for element in layer:
                if not element: continue
                print(f"Layer num: {element.layer_num}, index: {element.index}")
                print(f"phi: {np.round(element.phi, 4)}")

    def relax_mesh(self): # Absorb phase_diffs
        for phase in self.phantom_phases:
            self.print_elements()
            print(f"\n\ndissolving phase_shifter({phase.layer_num}, {phase.index}) with offset {phase.phi}")
            index = phase.index
            layer_num = phase.layer_num
            phi = -phase.phi
            for i in range(index+1, len(self.mesh_elements[0])):
                print(f"i = {i}")
                if self.mesh_elements[layer_num-1][i]:
                    self.mesh_elements[layer_num-1][i].add_phase_offset(phi)
                    print(f"   offset of {phi} added at [{layer_num-1}][{i}]")
                if self.mesh_elements[layer_num][i]:
                    self.mesh_elements[layer_num][i].add_phase_offset(-phi)
                    print(f"   offset of {-phi} added at [{layer_num}][{i}]")

    def is_unitary(self, m):
        m = np.matrix(m)
        return np.allclose(np.eye(m.shape[0]), m.H * m)


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
            i += 1
            if not element: continue
            if isinstance(element, Mzi): print(f"mzi({element.layer_num}, {element.index}): phi1 = {element.phi[0]}, phi2 = {element.phi[1]}")
            if isinstance(element, Phase_shifter): print(f"phase_shifter({element.layer_num}, {element.index}): phi = {element.phi}")


if __name__=="__main__":
    main()
