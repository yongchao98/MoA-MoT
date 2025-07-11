import math

def combinations(n, k):
    """Calculates the number of combinations C(n, k)."""
    if k < 0 or k > n:
        return 0
    # Use integer division for a clean result
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

# Part (b) parameters
k = 2
m = 5
t = 1

# A theorem for cross-intersecting multiset families states that for m >= k + t,
# the sum |F| + |G| is maximized by taking F and G to be the family of all
# k-multisets containing a fixed element. We calculate the size of this family.
# The size is the number of ways to choose the other k-1 elements from [m]
# with repetition, which is given by C(m + (k-1) - 1, k-1).

n_comb = m + k - 2
k_comb = k - 1
size_of_family = combinations(n_comb, k_comb)

# The sum-maximal value is twice the size of this family.
max_sum = 2 * size_of_family

print("Calculation for part (b):")
print(f"Given m = {m}, k = {k}, and t = {t}.")
print(f"The size of the family of all k-multisets containing a fixed element is C(m + k - 2, k - 1).")
print(f"This is C({m} + {k} - 2, {k} - 1) = C({n_comb}, {k_comb}).")
print(f"C({n_comb}, {k_comb}) = {size_of_family}.")
print(f"The maximum sum |F| + |G| is 2 * |A_i|, where A_i is such a family.")
print(f"Therefore, the maximal sum is 2 * {size_of_family} = {max_sum}.")
