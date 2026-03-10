import cmath as cm
import numpy as np

# 2x2 Mach-Zehnder Interferometer
phi = [cm.pi/2, 4*cm.pi/5]

s_in = np.array([1/cm.sqrt(2), 1/cm.sqrt(2)]).reshape((-1, 1))

sigma = (phi[0] + phi[1])/2
delta = (phi[0] - phi[1])/2

M = np.identity(2, dtype='complex')

m = np.array([[cm.exp(sigma*complex(0, 1)) * cm.sin(delta), cm.exp(sigma*complex(0,1)) * cm.cos(delta)], [cm.exp(sigma*complex(0, 1)) * cm.cos(delta), -cm.exp(sigma*complex(0, 1)) * cm.sin(delta)]])

x_offset = 0
y_offset = 0

M[x_offset:x_offset+m.shape[0], y_offset:y_offset+m.shape[1]] = m

print(M)
print(f"Input state: {s_in}")

res = np.matmul(M, s_in)
print(f"Resulting state:\n {res}")
