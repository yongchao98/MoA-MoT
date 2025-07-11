import numpy as np

# The initial K-matrix for the bosonic nu=2 state is the Pauli sigma_x matrix.
K_boson = np.array([[0, 1],
                    [1, 0]])

# The number of flux quanta attached to each fermion is 2.
# This corresponds to adding m * I to the K-matrix, where m=2 and I is the identity matrix.
m = 2
I = np.identity(2)
K_flux = m * I

# The K-matrix of the resulting fractional state is the sum of the bosonic K-matrix
# and the flux attachment term.
K_final = K_boson + K_flux

# Print the final result step by step
print("The K-matrix for the final electronic state is calculated as follows:")
print("K_electron = K_boson + m * I")
print(f"K_electron = [[{K_boson[0,0]}, {K_boson[0,1]}], [{K_boson[1,0]}, {K_boson[1,1]}]] + {m} * [[{int(I[0,0])}, {int(I[0,1])}], [{int(I[1,0])}, {int(I[1,1])}]]")
print(f"K_electron = [[{K_boson[0,0]}, {K_boson[0,1]}], [{K_boson[1,0]}, {K_boson[1,1]}]] + [[{int(K_flux[0,0])}, {int(K_flux[0,1])}], [{int(K_flux[1,0])}, {int(K_flux[1,1])}]]")
print("The resulting K-matrix is:")
print(f"[[{int(K_final[0,0])}, {int(K_final[0,1])}], [{int(K_final[1,0])}, {int(K_final[1,1])}]]")
