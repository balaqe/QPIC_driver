import numpy as np

class Clements:
    def __init__(self, n_ch: int):
        num_layers = n_ch

        self.M_layer = np.zeros((num_layers, n_ch, n_ch), dtype='complex')
        for layer in range(self.M_layer.shape[0]):
            self.M_layer[layer] = np.identity(n_ch, dtype='complex')

        self.eff_unitary = np.identity(n_ch, dtype='complex')
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

                m = np.array([[np.exp(self.sigma[layer_num][mzi]*1j) * np.sin(self.delta[layer_num][mzi]),\
                               np.exp(self.sigma[layer_num][mzi]*1j) * np.cos(self.delta[layer_num][mzi])],\
                              [np.exp(self.sigma[layer_num][mzi]*1j) * np.cos(self.delta[layer_num][mzi]),\
                               -np.exp(self.sigma[layer_num][mzi]*1j) * np.sin(self.delta[layer_num][mzi])]])

                self.M_layer[layer_num][x_offset:x_offset+m.shape[0], y_offset:y_offset+m.shape[1]] = m

            self.eff_unitary = np.matmul(self.M_layer[layer_num], self.eff_unitary)
    
    def execute(self, state: list[complex]):
        self.state = np.array(state, dtype='complex').reshape(-1, 1)
        self.state = np.matmul(self.eff_unitary, self.state)
        return self.state

    def display(self):
        print(self.eff_unitary)

chip = Clements(6)

chip.configure(layer_num=0, mzi=0, phi=[np.pi/2, 3*np.pi/2])
chip.configure(layer_num=0, mzi=1, phi=[0, 0])
chip.configure(layer_num=0, mzi=2, phi=[0, 0])

chip.configure(layer_num=1, mzi=0, phi=[np.pi/2, 3*np.pi/2])
chip.configure(layer_num=1, mzi=1, phi=[0, 0])

chip.configure(layer_num=2, mzi=0, phi=[np.pi/2, 3*np.pi/2])
chip.configure(layer_num=2, mzi=1, phi=[0, 0])
chip.configure(layer_num=2, mzi=2, phi=[0, 0])

chip.configure(layer_num=3, mzi=0, phi=[0, 0])
chip.configure(layer_num=3, mzi=1, phi=[0, 0])

chip.configure(layer_num=4, mzi=0, phi=[0, 0])
chip.configure(layer_num=4, mzi=1, phi=[0, 0])
chip.configure(layer_num=4, mzi=2, phi=[0, 0])

chip.configure(layer_num=5, mzi=0, phi=[np.pi/2, 3*np.pi/2])
chip.configure(layer_num=5, mzi=1, phi=[0, 0])
chip.build()

state = chip.execute([1, 0, 0, 0, 0, 0])
print("Resulting state:\n")
print(np.round(state, 2))