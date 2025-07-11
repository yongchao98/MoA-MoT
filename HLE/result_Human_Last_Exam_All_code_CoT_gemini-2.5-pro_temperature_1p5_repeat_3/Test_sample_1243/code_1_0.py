import math
from scipy.special import comb

# Part (a): Number of 2-isotropic vectors for Z^16
# This is the number of non-zero binary vectors of length 16 with weight divisible by 4.
n_a = 16
weights_a = [4, 8, 12, 16]
coeffs_a = [comb(n_a, k, exact=True) for k in weights_a]
ans_a = sum(coeffs_a)

# Part (b): Possibility of an even 3-neighbor for Z^8
# The reasoning leads to a contradiction, so the answer is no.
ans_b = "no"

# Part (c): Number of unimodular 2-neighbors of Z^12
# This requires x.x = 4. The vectors must have 4 non-zero entries (+/-1).
# For each set of 4 positions, there are 2 distinct lattices.
n_c = 12
k_c = 4
C12_4 = comb(n_c, k_c, exact=True)
ans_c = 2 * C12_4

# Print the final formatted answer as requested.
# For (a), the equation is the sum of binomial coefficients for weights 4, 8, 12, and 16.
print("Final Answer:")
print(f"(a) The number of vectors is C(16,4) + C(16,8) + C(16,12) + C(16,16) = {coeffs_a[0]} + {coeffs_a[1]} + {coeffs_a[2]} + {coeffs_a[3]} = {ans_a}")
# For (b), the answer is a simple 'no'.
print(f"(b) The answer is {ans_b}.")
# For (c), the equation is 2 times a binomial coefficient.
print(f"(c) The number of unimodular neighbors is 2 * C(12,4) = 2 * {C12_4} = {ans_c}")
