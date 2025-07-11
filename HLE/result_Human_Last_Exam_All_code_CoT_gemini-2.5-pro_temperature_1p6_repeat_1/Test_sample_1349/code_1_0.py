import math

# The supremum of X is given by the expression 44 / (36 * pi^2 + 25).
# Here, we calculate its numerical value.

pi = math.pi
numerator = 44
denominator_term1_coeff = 36
denominator_term2 = 25

# The minimum value of the Rayleigh quotient R is 1.8 * pi^2
R_min = 1.8 * pi**2
# The supremum of X is calculated using the formula X = 44 / (20*R + 25)
sup_X = numerator / (20 * R_min + denominator_term2)

# We can also substitute R_min = 1.8 * pi**2 directly
# X = 44 / (20 * 1.8 * pi**2 + 25) = 44 / (36 * pi**2 + 25)
calculated_denominator = denominator_term1_coeff * pi**2 + denominator_term2
sup_X_direct = numerator / calculated_denominator

print(f"The expression for the supremum of X is: {numerator} / ({denominator_term1_coeff} * pi^2 + {denominator_term2})")
print(f"The numerical value is: {sup_X_direct}")
