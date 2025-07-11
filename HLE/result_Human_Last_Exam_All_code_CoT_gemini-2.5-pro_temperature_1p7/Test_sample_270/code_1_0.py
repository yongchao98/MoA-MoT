import numpy as np

# Step 1: Define the matrices for the right-handed Dehn twists D_a and D_b.
# D_a is a twist along the curve 'a' (homology class [1, 0]).
# D_b is a twist along the curve 'b' (homology class [0, 1]).
# These are standard representations in SL(2, Z).
D_a_matrix = np.array([[1, 1], 
                       [0, 1]])

D_b_matrix = np.array([[1, 0], 
                       [-1, 1]])

# Step 2: Compute the matrix for the composite mapping class M = D_a * D_b.
M = np.dot(D_a_matrix, D_b_matrix)

# Step 3: Decompose the element (D_a*D_b)^9 into (D_a*D_b)^6 and (D_a*D_b)^3.
# We calculate the matrix representations for each part.
M_pow3 = np.linalg.matrix_power(M, 3)
M_pow6 = np.linalg.matrix_power(M, 6)

# The matrix for (D_a*D_b)^6 is the identity matrix.
# In the mapping class group, (D_a*D_b)^6 is equal to one inverse Dehn twist 
# around the boundary, D_d^-1.
# The fractional Dehn twist coefficient (FDTC) of D_d^k is k.
fdtc_of_pow6 = -1

# The matrix for (D_a*D_b)^3 is -I. This corresponds to the hyperelliptic involution,
# which is a pi-rotation of the torus.
# Its FDTC is the rotation angle (pi) divided by a full twist (2*pi).
fdtc_of_pow3 = 0.5 # which is pi / (2*pi)

# Step 4: The total FDTC is the sum of the coefficients of the parts.
# c((D_a*D_b)^9) = c((D_a*D_b)^6) + c((D_a*D_b)^3)
total_fdtc = fdtc_of_pow6 + fdtc_of_pow3

# Step 5: Print the final calculation step-by-step.
print("The fractional Dehn twist coefficient of (D_a o D_b)^9 is calculated by decomposing the element.")
print("We use the relation (D_a o D_b)^9 = (D_a o D_b)^6 o (D_a o D_b)^3.")
print("\nThe coefficient is the sum of the coefficients of the components:")
print(f"c((D_a o D_b)^9) = c((D_a o D_b)^6) + c((D_a o D_b)^3)")
print("\nFrom the theory of mapping class groups, we have:")
print(f"1. c((D_a o D_b)^6) is the coefficient of the boundary twist D_d^-1, which is {fdtc_of_pow6}.")
print(f"2. c((D_a o D_b)^3) is the coefficient of the hyperelliptic involution, which is {fdtc_of_pow3}.")
print("\nTherefore, the final equation for the coefficient is:")
print(f"{fdtc_of_pow6} + {fdtc_of_pow3} = {total_fdtc}")
