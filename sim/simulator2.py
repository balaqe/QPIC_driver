import numpy as np
from decomp.tests.arbitrary_q import *

class Simulator:
    def simulate(self, state, mesh_elements):
        self.dim = len(state)
        res = np.array(state).reshape(-1, 1)
        res_mat = np.identity(self.dim, dtype=np.complex128)
        for layer in reversed(mesh_elements):
            for element in layer:
                if not element: continue
                if isinstance(element, Phase_shifter):
                    mat = self.phase_op(element)
                    print(f"phase (layer {element.layer_num} index {element.index}")
                    print(np.round(mat, 2))
                    res = mat @ res
                    res_mat = res_mat @ mat
                elif isinstance(element, Mzi):
                    mat = self.mzi_op(element)
                    print(f"mzi (layer {element.layer_num} index {element.index}")
                    print(np.round(mat, 2))
                    res = mat @ res
                    res_mat = res_mat @ mat
                else:
                    print(f"ERROR!! Element type {type(element)} not supported")
        return res_mat

    def phase_op(self, shifter: Phase_shifter):
        res = np.identity(self.dim, dtype=np.complex128)
        index = shifter.index
        # print(f"shifter.phi = {shifter.phi} at layer {shifter.layer_num} and index {shifter.index}")
        # print(f"type(shifter): {type(shifter)}")
        res[index][index] = np.exp(1j*shifter.phi)
        return res

    def mzi_op(self, mzi: Mzi):
        res = np.identity(self.dim, dtype=np.complex128)
        index = mzi.index
        sigma = (mzi.phi[0] + mzi.phi[1]) / 2
        delta = (mzi.phi[0] - mzi.phi[1]) / 2

        res[index][index] = np.exp(1j*sigma) * np.sin(delta)
        res[index][index+1] = np.exp(1j*sigma) * np.cos(delta)
        res[index+1][index] = np.exp(1j*sigma) * np.cos(delta)
        res[index+1][index+1] = -np.exp(1j*sigma) * np.sin(delta)

        return res

