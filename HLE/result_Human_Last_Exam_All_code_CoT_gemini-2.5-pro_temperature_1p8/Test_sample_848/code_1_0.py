import math

# This script calculates the integer part of 10^4 * L, where L is the limit described in the problem.
# The analysis of the Diophantine equations reveals that solutions are grouped into two families,
# corresponding to two possible integer invariants, k=13 and k=5.

# The asymptotic density of solutions for each family is related to the largest root of a characteristic quadratic equation.

# For the family of solutions where k=13, the characteristic equation is x^2 - 13*x + 1 = 0.
# The largest root determines the growth rate of the solutions in this family.
r_13_num = 13
r_13_sqrt_term = 165  # This comes from 13^2 - 4
r_13 = (r_13_num + math.sqrt(r_13_sqrt_term)) / 2

# For the family of solutions where k=5, the characteristic equation is x^2 - 5*x + 1 = 0.
# The largest root determines the growth rate for this second family.
r_5_num = 5
r_5_sqrt_term = 21  # This comes from 5^2 - 4
r_5 = (r_5_num + math.sqrt(r_5_sqrt_term)) / 2

# The limit L = lim_{N->inf} F(N)/ln(N) is the sum of contributions from both families.
# Each family contributes 2/ln(r) to the limit, where r is its growth rate. The factor of 2
# accounts for the symmetric pairs (a,b) and (b,a).
# The natural logarithm ln(r) is equivalent to acosh(k/2) where k is the invariant.
limit_contrib_13 = 2 / math.log(r_13)
limit_contrib_5 = 2 / math.log(r_5)

L = limit_contrib_13 + limit_contrib_5

# The problem asks for the integer part of 10^4 * L.
result_factor = 10000
final_value = result_factor * L
integer_part = int(final_value)

# As requested, here is the final equation with all its numerical components.
print("The limit L is given by the formula:")
print(f"L = 2 / ln(({r_13_num} + sqrt({r_13_sqrt_term})) / 2) + 2 / ln(({r_5_num} + sqrt({r_5_sqrt_term})) / 2)")
print(f"\nThis evaluates to:")
print(f"L = {L}")
print(f"The quantity to compute is {result_factor} * L = {final_value}")

# Finally, we print the integer part of the result.
print(f"\nThe integer part of the final result is {integer_part}")