import numpy as np

# The K-matrix for the nu=2 Bosonic Integer Quantum Hall state.
K_boson = np.array([[0, 1],
                    [1, 0]])

# The correction matrix accounting for the internal structure of the bosons.
# The diagonal '1's account for the fermionic statistics of the underlying electrons.
# The off-diagonal '2's account for the mutual statistics from two attached fluxes.
delta_K = np.array([[1, 2],
                      [2, 1]])

# The K-matrix of the resulting fractional state of electrons is found by
# adding the correction matrix to the boson K-matrix.
K_electron = K_boson + delta_K

# Print the explanation and the final equation.
print("The K-matrix of the original BIQH state is:")
print(f"K_boson = \n{K_boson}\n")
print("The correction matrix accounting for the fermion statistics and attached flux is:")
print(f"delta_K = \n{delta_K}\n")
print("The K-matrix of the resulting fractional state of electrons is K_electron = K_boson + delta_K:")
print(f"K_electron = \n{K_boson[0,0]} {K_boson[0,1]}   {delta_K[0,0]} {delta_K[0,1]}   {K_electron[0,0]} {K_electron[0,1]}")
print(f"           {K_boson[1,0]} {K_boson[1,1]} + {delta_K[1,0]} {delta_K[1,1]} = {K_electron[1,0]} {K_electron[1,1]}")
