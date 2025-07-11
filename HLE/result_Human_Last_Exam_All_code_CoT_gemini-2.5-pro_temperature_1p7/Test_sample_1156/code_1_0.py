import math

# The problem is to find the normalised density of the invariant measure for the map
# T(x) = 1/sqrt(x) mod 1.

# The derived formula for the density is rho(x) = 1 / (2 * ln(2) * sqrt(x) * (1 + sqrt(x))).

# We will print this formula. The numbers in the final equation are 1 and 2,
# and the constant ln(2).

coeff_val = 2 * math.log(2)

print("The normalised density of the invariant measure is given by the formula:")
print("rho(x) = 1 / ( (2 * ln(2)) * sqrt(x) * (1 + sqrt(x)) )")
print("\nWhere:")
print(f"The number 2 is an integer constant.")
print(f"The number 1 is an integer constant.")
print(f"ln(2) is the natural logarithm of 2, approximately {math.log(2):.6f}.")
print(f"The full coefficient in the denominator, 2 * ln(2), is approximately {coeff_val:.6f}.")
