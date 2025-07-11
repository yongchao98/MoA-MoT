import math

# The invariant density is of the form rho(x) = C * x^p
# From our derivation, we found the values for C and p.
C = 0.5
p = -0.5

# We print the final equation for the normalised density.
# The formula is rho(x) = 1 / (2 * sqrt(x))
print("The normalised density of the invariant measure for the map T(x) = (1-sqrt(x))^2 is:")
print(f"rho(x) = {C} / x^{abs(p)}")
print("This can also be written as:")
print(f"rho(x) = 1 / (2 * sqrt(x))")

# Outputting the numbers in the final equation as requested.
denominator_coefficient = 2
print(f"\nThe numbers in the final simplified equation rho(x) = 1 / ({denominator_coefficient} * sqrt(x)) are:")
print(f"Numerator: 1")
print(f"Denominator coefficient: {denominator_coefficient}")