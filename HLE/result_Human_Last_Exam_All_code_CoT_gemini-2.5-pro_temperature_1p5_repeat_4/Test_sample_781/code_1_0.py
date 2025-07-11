import math

# The problem is to find the largest number n for a decomposition of a special
# continuum X, which has 5 distinguished points {a, b, c, d, e}.
num_points = 5

# Our analysis shows that a valid decomposition can be constructed by taking all
# subcontinua C(p, q) formed by pairs of these special points.
# The number of such pairs is the number of combinations of choosing 2 items
# from a set of 5, denoted C(5, 2).
k = 2

# We calculate C(n, k) = n! / (k! * (n-k)!)
n = num_points
c_nk = math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

# The final output prints each number in the equation, as requested.
n_fact = math.factorial(n)
k_fact = math.factorial(k)
n_minus_k_fact = math.factorial(n - k)

print(f"The problem reduces to finding the number of pairs of points from the set of {n} special points.")
print(f"This is a combinatorial calculation: C({n}, {k})")
print(f"C({n}, {k}) = {n}! / ({k}! * ({n}-{k})!)")
print(f"C({n}, {k}) = {n_fact} / ({k_fact} * {n_minus_k_fact})")
print(f"C({n}, {k}) = {c_nk}")
print(f"Thus, the largest number n is {c_nk}.")