import math

# Define the number of variables for the Boolean function
n = 4

# The problem is to find a(4), the maximal number of prime implicants
# for a Boolean function of 4 variables.
# The formula to calculate this is a(n) = C(n, k) * (2^k),
# where k = floor((2*n + 2) / 3).

# First, we calculate the value of k for n = 4.
k_numerator = 2 * n + 2
k_denominator = 3
k_float = k_numerator / k_denominator
k = math.floor(k_float)

# Next, we calculate the components for the formula a(4).
# Binomial coefficient C(n, k) or "n choose k"
combinations = math.comb(n, k)
# Power of 2
power_of_2 = 2**k

# Finally, we calculate the result a(4)
result = combinations * power_of_2

# Print the calculation steps
print(f"The formula for the maximal number of prime implicants a(n) is:")
print(f"a(n) = C(n, k) * (2^k), with k = floor((2n + 2) / 3)")
print(f"\nFor n = {n}:")
print(f"k = floor((2 * {n} + 2) / {k_denominator}) = floor({k_numerator} / {k_denominator}) = floor({k_float:.2f}) = {k}")
print(f"\nPlugging n = {n} and k = {k} into the formula:")
# The user wants each number in the final equation to be printed.
print(f"a({n}) = C({n}, {k}) * (2^{k})")
print(f"a({n}) = {combinations} * {power_of_2}")
print(f"a({n}) = {result}")