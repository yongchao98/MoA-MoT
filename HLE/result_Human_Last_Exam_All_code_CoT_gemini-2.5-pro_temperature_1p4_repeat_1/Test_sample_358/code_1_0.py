import numpy as np

# Step 1: Define the Cartan matrices for the principal blocks.
# For A_5 at p=2, from the Modular Atlas Project.
C_A5 = np.array([
    [4, 2, 2],
    [2, 2, 1],
    [2, 1, 2]
])

# For C_2 at p=2.
C_C2 = np.array([
    [2]
])

# Step 2: Calculate the sum of entries for each matrix.
sum_A5 = np.sum(C_A5)
sum_C2 = np.sum(C_C2)

# Step 3: Calculate the total sum for the product group G = A_5 x C_2.
# The sum of entries of the Kronecker product is the product of the sums.
total_sum = sum_A5 * sum_C2

# Step 4: Print the final equation with all its components as requested.
print(f"The Cartan matrix for the principal block of A5 in characteristic 2 is:\n{C_A5}")
print(f"The sum of its entries is: {sum_A5}")
print("\n")
print(f"The Cartan matrix for the principal block of C2 in characteristic 2 is:\n{C_C2}")
print(f"The sum of its entries is: {sum_C2}")
print("\n")
print("The sum of all entries in the Cartan matrix for G = A5 x C2 is the product of these sums.")
print(f"Total Sum = {sum_A5} * {sum_C2} = {total_sum}")

# The user might be interested in the full Cartan matrix for G as well.
C_G = np.kron(C_A5, C_C2)
print("\nThe full Cartan matrix for G = A5 x C2 is:")
print(C_G)
print(f"The sum of entries of this matrix is {np.sum(C_G)}, which confirms our result.")