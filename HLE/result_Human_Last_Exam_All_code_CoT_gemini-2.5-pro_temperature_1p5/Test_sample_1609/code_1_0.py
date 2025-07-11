import math

def combinations(n, k):
    """
    Helper function to calculate the binomial coefficient C(n, k),
    which is n! / (k! * (n-k)!).
    """
    if k < 0 or k > n:
        return 0
    # This calculation is efficient for small numbers like the ones in this problem.
    # It calculates (n * (n-1) * ... * (n-k+1)) / (k * (k-1) * ... * 1)
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

# The problem asks for a(4), where a(n) is the maximal number of prime implicants
# of a Boolean function of n variables.
n_val = 4

# For n <= 5, a(n) is given by the formula a(n) = C(2n, n), where C is the binomial coefficient "n choose k".
# We will calculate a(4) using this identity.
N = 2 * n_val
K = n_val

# Calculate the final result.
result = combinations(N, K)

# To fulfill the request, we will print out the equation with all numbers.
numerator_terms = list(range(N, N - K, -1))
denominator_terms = list(range(K, 0, -1))

print(f"The maximal number of prime implicants, a(n), for a Boolean function of n variables follows the identity a(n) = C(2n, n) for n <= 5.")
print(f"For n = {n_val}, the calculation is a(4) = C(2*{n_val}, {n_val}) = C({N}, {K}).")
print(f"The final equation and its solution are:")

# Format the equation string.
equation_str = f"a(4) = C({N}, {K}) = ({' * '.join(map(str, numerator_terms))}) / ({' * '.join(map(str, denominator_terms))}) = {result}"
print(equation_str)