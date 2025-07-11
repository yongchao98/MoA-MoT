# The exponent of N in the threshold for the sum's magnitude
exponent_threshold = 3/8

# The order of the moment used in Chebyshev's inequality
p = 4

# Based on harmonic analysis results and confirmed by a specific construction,
# the 4th moment of the maximal function scales with N.
# integral(M(x)^4 dx) is O(N^power_integral)
power_integral = 1

# From Chebyshev's inequality, we have:
# |X| <= N^(-p * exponent_threshold) * integral(M(x)^4 dx)
# |X| <= N^(-p * exponent_threshold) * C * N^power_integral
# So the exponent alpha is given by:
alpha = power_integral - p * exponent_threshold

term1 = power_integral
term2 = p * exponent_threshold

print(f"The calculation for alpha is based on the formula: alpha = power_integral - p * exponent_threshold")
print(f"The exponent in the integral bound is: {term1}")
print(f"The exponent from the threshold condition is: {p} * {exponent_threshold} = {term2}")
print(f"So, the final value of alpha is: {term1} - {term2} = {alpha}")

<<< -0.5 >>>