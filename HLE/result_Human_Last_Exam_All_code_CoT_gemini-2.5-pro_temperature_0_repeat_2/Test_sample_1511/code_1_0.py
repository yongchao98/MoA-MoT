import math

def combinations(n, k):
    """Calculates the number of combinations C(n, k)"""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

# Parameters for part (b)
m = 5
k = 2
t = 1

# The answers to (a) and (c) are derived from theoretical principles.
answer_a = "No"
answer_c = "Yes"

# For part (b), we calculate the maximal sum using the formula from the
# Ahlswede-Khachatrian theorem for multisets, as m >= k + t.
# The formula is 2 * C(m + k - t - 1, k - t).
n_comb = m + k - t - 1
k_comb = k - t
comb_val = combinations(n_comb, k_comb)
max_sum = 2 * comb_val
answer_b = max_sum

# Print the calculation for part (b) as requested.
print(f"The calculation for part (b) is based on the formula 2 * C(m+k-t-1, k-t).")
print(f"With m={m}, k={k}, t={t}, this becomes: 2 * C({m}+{k}-{t}-1, {k}-{t}) = 2 * C({n_comb}, {k_comb}) = 2 * {comb_val} = {max_sum}.")

# Print the final answer in the specified format.
print(f"\n(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")