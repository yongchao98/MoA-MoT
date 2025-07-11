# This script calculates the other critical exponent mentioned in the problem.
# The critical exponents p_k for the cone in R^3 are related to k-linear
# interactions and are given by the formula p_k = 2*k / (k-1).

# We are given that one critical exponent is 4. This corresponds to k=2 (bilinear case).
# p_2 = 2 * 2 / (2 - 1) = 4.

# The other critical exponent corresponds to the next case, k=3 (trilinear case).
k = 3
numerator = 2 * k
denominator = k - 1
result = numerator / denominator

# The problem asks to output the numbers in the final equation.
# Here is the calculation for k=3:
print(f"2 * {k} / ({k} - 1) = {numerator} / {denominator} = {int(result)}")