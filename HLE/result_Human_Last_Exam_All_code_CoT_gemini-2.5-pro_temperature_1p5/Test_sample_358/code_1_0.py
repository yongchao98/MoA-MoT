import numpy as np

# This plan calculates the sum of entries of the Cartan matrix for the principal block of G = A5 x C2 at characteristic p=2.

# 1. The group algebra kG is k(A5 x C2) which is isomorphic to kA5 (tensor) kC2.
# 2. The principal block of kG is the tensor product of the principal blocks of kA5 and kC2.
#    A5 is simple, so its group algebra at p=2 is a single block (the principal one).
#    C2 is a 2-group, so its group algebra at p=2 is also a single block.
# 3. The Cartan matrix of the tensor product is the Kronecker product of the individual Cartan matrices.
# 4. The sum of all entries in a Kronecker product of matrices is the product of the sums of entries of each matrix.
# 5. We calculate the sum of entries for the Cartan matrix of A5 and C2 separately.
#    The Cartan matrix C is related to the decomposition matrix D by C = D^T * D.

# Step 1: Calculate the sum of entries for the Cartan matrix of A5.
# The decomposition matrix D_A5 for A5 at p=2 maps its 5 ordinary irreducible characters
# to its 4 modular irreducible characters.
D_A5 = np.array([
    [1, 0, 0, 0],  # Ordinary char 1 (dim 1) -> Modular char 1 (dim 1)
    [1, 1, 0, 0],  # Ordinary char 2 (dim 3) -> Mod 1 + Mod 2 (dim 1+2)
    [1, 0, 1, 0],  # Ordinary char 3 (dim 3) -> Mod 1 + Mod 3 (dim 1+2)
    [0, 1, 1, 0],  # Ordinary char 4 (dim 4) -> Mod 2 + Mod 3 (dim 2+2)
    [1, 0, 0, 1]   # Ordinary char 5 (dim 5) -> Mod 1 + Mod 4 (dim 1+4)
])

# Calculate the Cartan matrix for A5.
C_A5 = D_A5.T @ D_A5

# Calculate the sum of its entries.
sum_A5 = np.sum(C_A5)

# Step 2: Calculate the sum of entries for the Cartan matrix of C2.
# The decomposition matrix D_C2 for C2 at p=2 maps its 2 ordinary irreducible characters
# to its single modular irreducible character.
D_C2 = np.array([
    [1],  # Trivial character -> Trivial modular character
    [1]   # Sign character -> Trivial modular character
])

# Calculate the Cartan matrix for C2.
C_C2 = D_C2.T @ D_C2

# Calculate the sum of its entries.
sum_C2 = np.sum(C_C2)

# Step 3: Compute the final result.
total_sum = sum_A5 * sum_C2

print("This script calculates the sum of all entries in the Cartan matrix for the principal block of G = A5 x C2 at characteristic 2.")
print("\nFirst, we find the sum of entries for the Cartan matrix of the principal block of A5.")
print(f"The sum for A5 is: {sum_A5}")
print("\nNext, we find the sum of entries for the Cartan matrix of the principal block of C2.")
print(f"The sum for C2 is: {sum_C2}")
print("\nThe total sum is the product of the individual sums.")
print(f"Final Equation: {sum_A5} * {sum_C2} = {total_sum}")
print(f"\nThe sum of all the entries in the Cartan matrix of the principal block of k(A5 x C2) is {total_sum}.")

<<<34>>>