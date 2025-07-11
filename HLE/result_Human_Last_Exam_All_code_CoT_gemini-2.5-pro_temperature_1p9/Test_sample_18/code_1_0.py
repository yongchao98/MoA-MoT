import numpy as np

# The initial K-matrix for the nu=2 Bosonic Integer Quantum Hall state (sigma_x)
K_b = np.array([[0, 1], 
                [1, 0]])

# The identity matrix
I = np.identity(2, dtype=int)

# The transformation involves attaching 2 fluxes and pairing (changing statistics),
# which corresponds to an effective addition of 3 statistical fluxes.
# K_f = K_b + 3*I
m = 3
K_f = K_b + m * I

# Print the final result, showing each number in the equation as requested.
print("The K-matrix of the resulting fractional state is K_f.")
print("It is calculated as K_f = K_b + 3*I, where K_b is the initial bosonic K-matrix.")
print("\nFinal K-matrix, K_f:")
print(K_f)
print("\n--- Derivation ---")
print("K_f = ")
print(K_b)
print(" + {} * ".format(m))
print(I)
print(" = ")
print(K_f)
