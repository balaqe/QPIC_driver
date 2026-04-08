import numpy as np

class Clements:
    def __init__(self, n_ch: int):
        num_layers = n_ch+1

        self.M_layer = np.zeros((num_layers, n_ch, n_ch), dtype='complex')
        for layer in range(self.M_layer.shape[0]):
            self.M_layer[layer] = np.identity(n_ch, dtype='complex')

        self.eff_unitary = np.identity(n_ch, dtype='complex')
        self.phi = np.zeros((num_layers, n_ch//2, 2), dtype='complex') # Layer, MZI per layer, heater per MZI
        self.sigma = np.zeros((num_layers, n_ch//2), dtype='complex')
        self.delta = np.zeros((num_layers, n_ch//2), dtype='complex')
        self.shifter = np.empty((num_layers, n_ch), dtype='complex')
        self.shifter[:] = np.nan

    def configure_mzi(self, layer_num: int, mzi: int, phi: list[complex]):
        self.phi[layer_num][mzi] = np.array(phi)
        self.sigma[layer_num][mzi] = (phi[0] + phi[1]) / 2
        self.delta[layer_num][mzi] = (phi[0] - phi[1]) / 2

    def configure_phase(self, layer_num: int, shifter: int, phi: complex):
        self.shifter[layer_num][shifter] = np.exp(1j*phi)

    def build(self):
        for layer_num in range(self.sigma.shape[0]):
            for shifter in range(self.shifter.shape[1]):
                if not np.isnan(self.shifter[layer_num][shifter]):
                    self.M_layer[layer_num][shifter][shifter] = self.shifter[layer_num][shifter]

            if layer_num == 0: continue # Skip initial external shifter layer
        
            for mzi in range(self.sigma.shape[1]):
                x_offset = 0
                y_offset = 0
                if layer_num%2 == 0: # Calculate start offset
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
                # print(f"Layer:\n{self.M_layer[layer_num]}")
                # print(f"layer_num = {layer_num}; mzi = {mzi}")
                # print(f"Sigma: {self.sigma[layer_num][mzi]}, Delta: {self.delta[layer_num][mzi]}")

            self.eff_unitary = self.M_layer[layer_num] @ self.eff_unitary
    
    def execute(self, state: list[complex]):
        self.state = np.array(state, dtype='complex').reshape(-1, 1)
        self.state = np.matmul(self.eff_unitary, self.state)
        return self.state

    def display(self):
        print(self.eff_unitary)

chip = Clements(6)

chip.configure_phase(layer_num=0, shifter=0, phi=0)
chip.configure_phase(layer_num=0, shifter=1, phi=0)
chip.configure_phase(layer_num=0, shifter=2, phi=0)
chip.configure_phase(layer_num=0, shifter=3, phi=0)
chip.configure_phase(layer_num=0, shifter=4, phi=0)
chip.configure_phase(layer_num=0, shifter=5, phi=0)

chip.configure_mzi(layer_num=1, mzi=0, phi=[3*np.pi/2, np.pi/2])
chip.configure_mzi(layer_num=1, mzi=1, phi=[0, 0])
chip.configure_mzi(layer_num=1, mzi=2, phi=[0, 0])

chip.configure_mzi(layer_num=2, mzi=0, phi=[0, 0])
chip.configure_mzi(layer_num=2, mzi=1, phi=[0, 0])

chip.configure_mzi(layer_num=3, mzi=0, phi=[3*np.pi/2, np.pi/2])
chip.configure_mzi(layer_num=3, mzi=1, phi=[0, 0])
chip.configure_mzi(layer_num=3, mzi=2, phi=[0, 0])

chip.configure_mzi(layer_num=4, mzi=0, phi=[0, 0])
chip.configure_mzi(layer_num=4, mzi=1, phi=[0, 0])

chip.configure_mzi(layer_num=5, mzi=0, phi=[0, 0])
chip.configure_mzi(layer_num=5, mzi=1, phi=[0, 0])
chip.configure_mzi(layer_num=5, mzi=2, phi=[0, 0])

chip.configure_mzi(layer_num=6, mzi=0, phi=[3*np.pi/2, np.pi/2])
chip.configure_mzi(layer_num=6, mzi=1, phi=[0, 0])

# chip = Clements(2)

#chip.configure_mzi(layer_num=1, mzi=0, phi=[3*np.pi/2, np.pi/2])
# chip.configure_mzi(layer_num=1, mzi=0, phi=[0, 0])

chip.build()

state = chip.execute([1, 0, 0, 0, 0, 0])
# state = chip.execute([1, 0])

print("Resulting state:\n")
print(np.round(state, 2))
