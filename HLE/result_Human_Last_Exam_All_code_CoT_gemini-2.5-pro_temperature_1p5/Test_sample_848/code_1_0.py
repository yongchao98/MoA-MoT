import math

# This plan calculates the value based on the derived formula for the limit.
# The limit C is given by the sum of contributions from two families of solutions,
# corresponding to k=5 and k=13.
# C = 2/ln(alpha_5) + 2/ln(alpha_13)

# For the k=5 family, the characteristic equation is r^2 - 5r + 1 = 0.
# The largest root is alpha_5 = (5 + sqrt(5^2 - 4))/2 = (5 + sqrt(21))/2.
k5 = 5
alpha_5 = (k5 + math.sqrt(k5**2 - 4)) / 2

# For the k=13 family, the characteristic equation is r^2 - 13r + 1 = 0.
# The largest root is alpha_13 = (13 + sqrt(13^2 - 4))/2 = (13 + sqrt(165))/2.
k13 = 13
alpha_13 = (k13 + math.sqrt(k13**2 - 4)) / 2

# The limit C = lim_{N->inf} F(N)/ln(N)
limit_C = 2 / math.log(alpha_5) + 2 / math.log(alpha_13)

# The problem asks for the integer part of 10^4 * C.
final_value = 10000 * limit_C
integer_part = math.floor(final_value)

print(integer_part)
