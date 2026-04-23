import numpy as np
import sympy

u11 = sympy.Symbol('u11')
u12 = sympy.Symbol('u12')
u21 = sympy.Symbol('u21')
u22 = sympy.Symbol('u22')

U = np.array([[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]])

if U[0][3] != 0:
    sigma = np.angle(U[0][-1])
    # print(f"U[0][2] = {U[0][2]}; U[0][3] = {U[0][3]}")
    delta = sympy.atan(U[0][2]/U[0][3])
    phase_diff = np.angle(U[0][-1]) - np.angle(U[-1][0])

    print(f"sigma = {sigma}; delta = {delta}; phase_diff = {phase_diff}")

    c = []
    for i in range(4):
        c.append(sympy.Symbol(f'c{i}'))

    m = np.array([[sympy.exp(sigma*1j) * sympy.sin(delta),\
                    sympy.exp(sigma*1j) * sympy.cos(delta)],\
                    [sympy.exp(sigma*1j) * sympy.cos(delta),\
                    -sympy.exp(sigma*1j) * sympy.sin(delta)]])

    M = np.identity(4)
    x_offset = 2
    y_offset = 2
    M[x_offset:x_offset+m.shape[0], y_offset:y_offset+m.shape[1]] = m

    M[-1][0] *= sympy.exp(phase_diff*1j)

    print(f"M:\n{M}\n")
    print(f"j = {1j}, -j = {-1j}")

    c = np.array(c).reshape(-1, 1)

    res1 = U @ c
    print("\nresult1:")
    print(res1)

    res2 = M @ c
    print("\nresult2:")
    print(res2)
