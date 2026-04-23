import numpy as np
import sympy

u11 = sympy.Symbol('u11')
u12 = sympy.Symbol('u12')
u21 = sympy.Symbol('u21')
u22 = sympy.Symbol('u22')

# U = np.array([[u11, u12], [u21, u22]])
# U = 1/np.sqrt(2) * np.array([[1, 1], [1, -1]])
# U = np.array([[0, 1], [1, 0]])
U = np.array([[0, -1j], [1j, 0]])

# sigma = sympy.Symbol('sigma')
# delta = sympy.Symbol('delta')
sigma = np.angle(U[0][1])
delta = sympy.atan(U[0][0]/U[0][1])
phase_diff = np.angle(U[0][1]) - np.angle(U[1][0])

c = []

for i in range(2):
    c.append(sympy.Symbol(f'c{i}'))

M = np.array([[sympy.exp(sigma*1j) * sympy.sin(delta),\
                sympy.exp(sigma*1j) * sympy.cos(delta)],\
                [sympy.exp(sigma*1j) * sympy.cos(delta),\
                -sympy.exp(sigma*1j) * sympy.sin(delta)]])

M = np.array([[1, 0], [0, sympy.exp(phase_diff*1j)]]) @ M # Add phase difference


print(f"M:\n{M}\n")
print(f"j = {1j}, -j = {-1j}")

c = np.array(c).reshape(-1, 1)

res1 = U @ c
print("\nresult1:")
print(res1)

# c[0] *= sympy.exp(phase_diff*1j)

res2 = M @ c
print("\nresult2:")
print(res2)
