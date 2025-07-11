import math

def combinations(n, k):
    """Computes the binomial coefficient C(n, k)."""
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

# --- Part (a) ---
# For L = Z^16, we count the number of non-zero vectors x_bar in (Z/2Z)^16
# such that the weight (number of 1s) is a multiple of 4.
# The possible weights are 4, 8, 12, 16.
n_a = 16
c16_4 = combinations(n_a, 4)
c16_8 = combinations(n_a, 8)
c16_12 = combinations(n_a, 12)
c16_16 = combinations(n_a, 16)
answer_a = c16_4 + c16_8 + c16_12 + c16_16

# --- Part (b) ---
# As explained in the thinking steps, the conditions on the vector x lead to a contradiction.
# Specifically, x must be a constant vector c*(1,...,1) where c is not a multiple of 3.
# However, the condition for the neighbor lattice to be even implies that c must be a
# multiple of 3. Thus, it's not possible.
answer_b = "no"

# --- Part (c) ---
# The number of unimodular 2-neighbors of Z^12 is 2 times the number of
# non-zero vectors x_bar in (Z/2Z)^12 with weight divisible by 4.
# The possible weights are 4, 8, 12.
n_c = 12
c12_4 = combinations(n_c, 4)
c12_8 = combinations(n_c, 8)
c12_12 = combinations(n_c, 12)
num_vectors_c = c12_4 + c12_8 + c12_12
answer_c = 2 * num_vectors_c

# --- Final Output ---
print(f"(a) {c16_4} + {c16_8} + {c16_12} + {c16_16} = {answer_a}")
print(f"(b) {answer_b}")
print(f"(c) 2 * ({c12_4} + {c12_8} + {c12_12}) = {answer_c}")

# The final answer should be provided in the required format.
# print(f"({answer_a}; {answer_b}; {answer_c})")
