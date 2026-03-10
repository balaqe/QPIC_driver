import cmath as cm
import numpy as np

# 2x2 Mach-Zehnder Interferometer
n_ch = 6
phi = [cm.pi/2, 4*cm.pi/5]

assert n_ch % 2 == 0


s_in = np.zeros(n_ch).reshape((-1, 1))
s_in[0] = 1

sigma = (phi[0] + phi[1])/2
delta = (phi[0] - phi[1])/2

M = []

M.append(np.identity(n_ch, dtype='complex'))

m = np.array([[cm.exp(sigma*complex(0, 1)) * cm.sin(delta), cm.exp(sigma*complex(0,1)) * cm.cos(delta)], [cm.exp(sigma*complex(0, 1)) * cm.cos(delta), -cm.exp(sigma*complex(0, 1)) * cm.sin(delta)]])

x_offset = 0
y_offset = 0

M[0][x_offset:x_offset+m.shape[0], y_offset:y_offset+m.shape[1]] = m

#print(M[0])
#print(f"Input state: {s_in}")

res = np.matmul(M[0], s_in)
#print(f"Resulting state:\n {res}")

class Clements:
    def __init__(self, n_ch: int):
        num_layers = n_ch
        self.M_layer = np.zeros((num_layers, n_ch, n_ch), dtype='complex')
        self.phi = np.zeros((n_ch, n_ch//2, 2), dtype='complex') # Layer, MZI per layer, heater per MZI
        self.sigma = np.zeros((n_ch, n_ch//2))
        self.delta = np.zeros((n_ch, n_ch//2))

    def configure(self, layer_num: int, mzi: int, phi: list[complex]):
        self.phi[layer_num][mzi] = np.array(phi)
        self.sigma[layer_num][mzi] = (phi[0] + phi[1]) / 2
        self.delta[layer_num][mzi] = (phi[0] - phi[1]) / 2

    def build(self):
        for layer_num in range(self.sigma.shape[0]):
            for mzi in range(self.sigma.shape[1]):
                x_offset = 0
                y_offset = 0
                if layer_num%2 != 0: # Calculate start offset
                    x_offset = 1
                    y_offset = 1
                    if mzi >= (self.sigma.shape[1] - 1):
                        continue

                x_offset += 2*mzi
                y_offset += 2*mzi

                m = np.array([[cm.exp(self.sigma[layer_num][mzi]*complex(0,1)) * cm.sin(self.delta[layer_num][mzi]),\
                            cm.exp(self.sigma[layer_num][mzi]*complex(0,1)) * cm.cos(self.delta[layer_num][mzi])],\
                            [cm.exp(self.sigma[layer_num][mzi]*complex(0, 1)) * cm.cos(self.delta[layer_num][mzi]),\
                            -cm.exp(self.sigma[layer_num][mzi]*complex(0, 1)) * cm.sin(self.delta[layer_num][mzi])]])

                self.M_layer[layer_num][x_offset:x_offset+m.shape[0], y_offset:y_offset+m.shape[1]] = m

    def display(self):
        print(self.M_layer)

chip = Clements(6)
chip.configure(layer_num=0, mzi=0, phi=[cm.pi, cm.pi/2])
chip.build()
chip.display()
