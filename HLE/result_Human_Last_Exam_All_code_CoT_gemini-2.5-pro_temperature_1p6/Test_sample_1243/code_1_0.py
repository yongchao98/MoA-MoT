import math

# (a) If n = 16 and L = Z^16, how many distinct 2-isotropic vectors can be found
# such that they satisfy the conditions for defining a 2-neighbor N of L?
# This means we count the number of non-zero vectors x_bar in (Z/2Z)^16
# such that the Hamming weight w(x_bar) is a multiple of 4.
# The total number of such vectors (including the zero vector) is S0 = 2^(n-2) + 2^(n/2-1)*cos(n*pi/4)
n_a = 16
# cos(16*pi/4) = cos(4*pi) = 1
term1_a = 2**(n_a - 2)
term2_a = 2**(n_a // 2 - 1)
ans_a = term1_a + term2_a - 1
print(f"For (a), the number is {term1_a} + {term2_a} - 1 = {ans_a}")

# (b) Let x in Z^8 be a 3-primitive vector. Is it possible for x to be 3-isotropic
# in such a way that the resulting 3-neighbor N_3(x) is even?
# As explained in the reasoning, this is not possible because the sublattice M_3(x)
# can never be an even lattice.
ans_b = "no"
print(f"For (b), the answer is {ans_b}")

# (c) Determine the number of unimodular 2-neighbors N of Z^12.
# This corresponds to counting non-zero vectors in (Z/2Z)^12 with weight divisible by 4.
n_c = 12
# cos(12*pi/4) = cos(3*pi) = -1
term1_c = 2**(n_c - 2)
term2_c = 2**(n_c // 2 - 1)
ans_c = term1_c - term2_c - 1
print(f"For (c), the number is {term1_c} - {term2_c} - 1 = {ans_c}")

# Format the final answer string for clarity, though the individual results are printed above.
final_answer = f"(a) {ans_a}; (b) {ans_b}; (c) {ans_c}"
print(f"\nFinal Answer: {final_answer}")