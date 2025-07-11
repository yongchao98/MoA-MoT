import math

# Part (a): Number of 2-isotropic vectors for a 2-neighbor of Z^16
n_a = 16
# These are non-zero vectors in (Z/2Z)^16 with weight divisible by 4.
# Possible weights are 4, 8, 12, 16.
c_a_4 = math.comb(n_a, 4)
c_a_8 = math.comb(n_a, 8)
c_a_12 = math.comb(n_a, 12)
c_a_16 = math.comb(n_a, 16)
total_a = c_a_4 + c_a_8 + c_a_12 + c_a_16

print("(a) The number of distinct 2-isotropic vectors is given by the sum:")
print(f"C(16, 4) + C(16, 8) + C(16, 12) + C(16, 16)")
print(f"which evaluates to {c_a_4} + {c_a_8} + {c_a_12} + {c_a_16} = {total_a}.")
print("-" * 20)

# Part (b): Possibility of an even 3-neighbor of Z^8
print("(b) The answer is based on a theoretical argument.")
print("It is not possible for the neighbor lattice to be even, as the required sublattice L_3(x) can never be even.")
print("The answer is no.")
print("-" * 20)

# Part (c): Number of unimodular 2-neighbors of Z^12
n_c = 12
k_c = 4
# The number of such neighbors corresponds to the number of binary vectors
# of length 12 with weight 4.
total_c = math.comb(n_c, k_c)
print("(c) The number of unimodular 2-neighbors is C(12, 4).")
print(f"C({n_c}, {k_c}) = {total_c}.")
