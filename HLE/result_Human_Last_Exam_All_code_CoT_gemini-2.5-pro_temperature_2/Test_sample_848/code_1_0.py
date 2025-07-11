import math

# The problem is to find the integer part of 10^4 * lim_{N->inf} F(N)/ln(N).
# The limit is given by the sum of contributions from two families of solutions,
# corresponding to k=13 and k=5.

# For k=13, the growth rate is determined by phi_13 = (13 + sqrt(165))/2
k1 = 13
phi1_numerator = k1 + math.sqrt(k1**2 - 4)
phi1 = phi1_numerator / 2
term1 = 2 / math.log(phi1)

# For k=5, the growth rate is determined by phi_5 = (5 + sqrt(21))/2
k2 = 5
phi2_numerator = k2 + math.sqrt(k2**2 - 4)
phi2 = phi2_numerator / 2
term2 = 2 / math.log(phi2)

# The limit is the sum of these two terms.
limit_val = term1 + term2

# We need to find the integer part of 10^4 times this limit.
result = 10000 * limit_val
integer_part = int(result)

print("The first part of the sum in the limit expression is 2 / ln((13 + sqrt(165))/2)")
print(f"Value: {term1}")
print("The second part of the sum in the limit expression is 2 / ln((5 + sqrt(21))/2)")
print(f"Value: {term2}")
print(f"The total limit value is the sum: {limit_val}")
print(f"The final result is 10^4 * limit: {result}")
print(f"The integer part of the final result is: {integer_part}")
