import math

# Calculate the normalization constant C = 1 / (2 * ln(2))
C = 1 / (2 * math.log(2))

# Print the formula for the normalized density function
# We use string formatting to include the calculated value of C
print("The normalised density of the invariant measure is rho(x) = C / (x + sqrt(x))")
print(f"where the constant C is 1/(2*ln(2))")
print("\nThe equation for the density is:")
print(f"rho(x) = {C:.8f} / (x + sqrt(x))")
