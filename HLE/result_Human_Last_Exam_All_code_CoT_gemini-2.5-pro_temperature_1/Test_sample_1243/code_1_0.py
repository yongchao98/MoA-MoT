import math

def combinations(n, k):
    """Calculates the binomial coefficient C(n, k)"""
    if k < 0 or k > n:
        return 0
    # Use integer division for precise results
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

# --- Part (a) ---
# Count 2-primitive, 2-isotropic vectors in Z^16. This corresponds to non-zero
# vectors in (Z/2Z)^16 whose weight k is a multiple of 4.
n_a = 16
k_values_a = [4, 8, 12, 16]
c16_4 = combinations(n_a, 4)
c16_8 = combinations(n_a, 8)
c16_12 = combinations(n_a, 12)
c16_16 = combinations(n_a, 16)
answer_a = c16_4 + c16_8 + c16_12 + c16_16
print("(a) The number of vectors is the sum of the ways to choose k odd components where k is a multiple of 4:")
print("C(16, 4) + C(16, 8) + C(16, 12) + C(16, 16) = {} + {} + {} + {} = {}".format(c16_4, c16_8, c16_12, c16_16, answer_a))


# --- Part (b) ---
# Determine if a 3-neighbor of Z^8 can be even.
answer_b = "no"
print("\n(b) The answer is no. Any such neighbor N must contain a vector with an odd norm (e.g., a standard basis vector e_j), so it cannot be an even lattice.")

# --- Part (c) ---
# Count the number of unimodular 2-neighbors of Z^12.
# This is twice the number of valid sublattices M_2(x).
n_c = 12
k_values_c = [4, 8, 12]
c12_4 = combinations(n_c, 4)
c12_8 = combinations(n_c, 8)
c12_12 = combinations(n_c, 12)
num_sublattices = c12_4 + c12_8 + c12_12
answer_c = num_sublattices * 2
print("\n(c) The number of neighbors is 2 times the number of valid sublattices M_2(x).")
print("Number of sublattices = C(12, 4) + C(12, 8) + C(12, 12) = {} + {} + {} = {}.".format(c12_4, c12_8, c12_12, num_sublattices))
print("Total number of neighbors = ({} + {} + {}) * 2 = {} * 2 = {}".format(c12_4, c12_8, c12_12, num_sublattices, answer_c))

print("\n---")
# Format the final answer as requested
final_answer = "(a) [{}]; (b) [{}]; (c) [{}]".format(answer_a, answer_b, answer_c)
print(final_answer)
print("<<<{}>>>".format(final_answer))