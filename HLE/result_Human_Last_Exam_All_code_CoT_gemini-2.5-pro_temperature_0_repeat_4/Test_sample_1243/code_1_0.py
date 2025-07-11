import math

def nCr(n, r):
    """Calculates the binomial coefficient 'n choose r'."""
    if r < 0 or r > n:
        return 0
    # Using math.comb for efficiency and clarity (available in Python 3.8+)
    # For older versions, one could use math.factorial:
    # f = math.factorial
    # return f(n) // f(r) // f(n-r)
    return math.comb(n, r)

# --- Part (a) ---
# The number of 2-isotropic, 2-primitive vectors for n=16 corresponds to
# the number of non-zero vectors in (Z/2Z)^16 whose Hamming weight is a multiple of 4.
n_a = 16
c16_4 = nCr(n_a, 4)
c16_8 = nCr(n_a, 8)
c16_12 = nCr(n_a, 12)
c16_16 = nCr(n_a, 16)
ans_a = c16_4 + c16_8 + c16_12 + c16_16

# --- Part (b) ---
# Based on the logical argument in the plan, it's impossible for such a neighbor to be even.
ans_b = "no"

# --- Part (c) ---
# The number of unimodular 2-neighbors of Z^12 corresponds to the number of
# non-zero vectors in (Z/2Z)^12 whose Hamming weight is a multiple of 4.
n_c = 12
c12_4 = nCr(n_c, 4)
c12_8 = nCr(n_c, 8)
c12_12 = nCr(n_c, 12)
ans_c = c12_4 + c12_8 + c12_12

# --- Final Output ---
# Print the results in the required format.
# The prompt asks to output each number in the final equation.

print("Calculations:")
print(f"(a) The number of vectors is C(16, 4) + C(16, 8) + C(16, 12) + C(16, 16) = {c16_4} + {c16_8} + {c16_12} + {c16_16} = {ans_a}")
print(f"(b) The answer is '{ans_b}' based on a combinatorial contradiction.")
print(f"(c) The number of neighbors is C(12, 4) + C(12, 8) + C(12, 12) = {c12_4} + {c12_8} + {c12_12} = {ans_c}")
print("\nFinal Answer:")
print(f"(a) {ans_a}; (b) {ans_b}; (c) {ans_c}")
print(f"<<<(a) {ans_a}; (b) {ans_b}; (c) {ans_c}>>>")