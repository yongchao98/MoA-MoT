import numpy as np

# Step 1: Determine the Cartan matrix for the principal block of kA_5 (p=2).
# This is equivalent to finding the Cartan matrix for kA_4 (p=2).
# The group A_4 has 4 ordinary irreducible characters and 3 simple modules in characteristic 2.
# The decomposition matrix D for kA_4 is known to be:
D_A4 = np.array([[1, 0, 0],
                 [0, 1, 0],
                 [0, 0, 1],
                 [1, 1, 1]])

# The Cartan matrix C is D^T * D.
C_A5_principal = D_A4.T @ D_A4

# Step 2: Determine the Cartan matrix for the principal block of kC_2 (p=2).
# kC_2 has one simple module, and its projective cover has length 2.
C_C2_principal = np.array([[2]])

# Step 3: Compute the Cartan matrix for G = A_5 x C_2.
# This is the Kronecker product of the individual Cartan matrices.
C_G = np.kron(C_A5_principal, C_C2_principal)

# Step 4: Calculate and print the sum of all entries.
total_sum = np.sum(C_G)

print("The Cartan matrix for the principal block of k(A_5 x C_2) in characteristic 2 is C_G, which is given by:")
print(C_G)
print("\nThe sum of all entries is found by adding every element in this matrix.")

# Format the final equation string as requested.
entries_str = [str(entry) for entry in C_G.flatten()]
equation_str = " + ".join(entries_str)

print(f"The sum is: {equation_str} = {int(total_sum)}")
