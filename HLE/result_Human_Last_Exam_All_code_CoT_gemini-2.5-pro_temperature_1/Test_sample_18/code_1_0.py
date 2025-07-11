import numpy as np

# The K-matrix for the nu=2 Bosonic Integer Quantum Hall state is the Pauli matrix sigma_x.
K_boson = np.array([[0, 1],
                    [1, 0]])

# The transformation rule between a fermionic K-matrix (K_fermion) and the
# K-matrix of its Cooper-paired bosonic state (K_boson) is K_boson = 4 * K_fermion.
# To find the K-matrix of the composite fermions, we reverse this transformation.
# K_fermion = K_boson / 4
K_fermion = K_boson / 4.0

# Print the final equation and the resulting K-matrix
print("The K-matrix of the bosonic state (Cooper pairs) is K_boson:")
print(K_boson)
print("\nThe transformation to the constituent fermion state is: K_fermion = K_boson / 4")
print("\nThe resulting K-matrix of the fractional state (the composite fermions) is:")
print(K_fermion)
