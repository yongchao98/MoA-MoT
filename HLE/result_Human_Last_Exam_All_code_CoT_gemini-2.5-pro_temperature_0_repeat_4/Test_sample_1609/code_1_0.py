import math

def combinations(n, k):
    """
    Helper function to calculate the binomial coefficient C(n, k),
    also known as "n choose k".
    """
    if k < 0 or k > n:
        return 0
    # math.comb is available in Python 3.8+
    return math.comb(n, k)

# The maximal number of prime implicants, a(n), is calculated using the recurrence:
# a(n) = 2 * a(n-1) + C(n-1, floor(n/2))
# with the base case a(0) = 0.

# We compute the sequence up to n=4.
a = {}
a[0] = 0

# Calculate a(1), a(2), a(3)
a[1] = 2 * a[0] + combinations(0, math.floor(1 / 2)) # 2*0 + C(0,0) = 1
a[2] = 2 * a[1] + combinations(1, math.floor(2 / 2)) # 2*1 + C(1,1) = 3
a[3] = 2 * a[2] + combinations(2, math.floor(3 / 2)) # 2*3 + C(2,1) = 8

# Now, calculate the final value for a(4)
n_final = 4
prev_a = a[n_final - 1]
n_comb = n_final - 1
k_comb = math.floor(n_final / 2)
comb_val = combinations(n_comb, k_comb)
a[n_final] = 2 * prev_a + comb_val

# Print the final equation showing each number in the calculation
print(f"a(4) = 2 * a(3) + C({n_comb}, {k_comb}) = 2 * {prev_a} + {comb_val} = {a[n_final]}")