import math

# The problem asks for a(4), the maximal number of prime implicants
# for a Boolean function of n=4 variables.
n = 4

# For n >= 3, the value is given by the formula:
# a(n) = n * C(n-1, floor((n-1)/2))
# where C is the binomial coefficient (combinations).

# Calculate the terms for the formula
n_minus_1 = n - 1
k = math.floor(n_minus_1 / 2)

# Calculate the binomial coefficient C(n-1, k)
# C(3, 1) in this case
combinations = math.comb(n_minus_1, k)

# Calculate the final result for a(4)
result = n * combinations

# Print the final equation with all the numbers
print(f"a({n}) = {n} * C({n_minus_1}, {k}) = {result}")
