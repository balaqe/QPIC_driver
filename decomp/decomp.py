import numpy as np
import math
import cmath as cm

def decomp(U):
    V = U # Set auxillary matrix
    n = np.shape(U)[0]
    m = np.shape(U)[1]
    print(f"n = {n}   m = {m}")
    for j in range(m-1): # row = 0 to m-2
        if j%2 == 0: # If iter is even (odd in paper due to one-based indexing)
            x = m - 1 # Offset due to zero-based indexing
            y = j
            P = np.identity(m, dtype='complex')
            phi = cm.phase(P[x, y]) - cm.phase(P[x, y+1])
            P[j, j] = np.exp(phi*1j)
            V = np.matmul(V, P)
            # print(f"x = {x}    y = {y}")
            # print(f"V = {V}")
            for k in range(j): # k = 0 to j-1
                M = np.identity(m, dtype='complex')
                delta = 0
                print(f"x={x}, y={y}, j={j}, k={k}")
                #V = np.matmul(V, M)

U = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
print(U)
decomp(U)
