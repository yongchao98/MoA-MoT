# The ratio of the maximum tangential stress (σ_t_max) to the internal
# pressure (P_i) in a thick-walled cylinder is given by the formula:
# Ratio = (r_o^2 + r_i^2) / (r_o^2 - r_i^2)
# where r_o is the outer radius and r_i is the inner radius.

# We are given that the outer radius is twice the inner radius,
# so the ratio of radii, k = r_o / r_i, is 2.
k = 2

# We can rewrite the formula in terms of k:
# Ratio = ( (k*r_i)^2 + r_i^2 ) / ( (k*r_i)^2 - r_i^2 )
# Ratio = ( k^2*r_i^2 + r_i^2 ) / ( k^2*r_i^2 - r_i^2 )
# Ratio = ( r_i^2 * (k^2 + 1) ) / ( r_i^2 * (k^2 - 1) )
# Ratio = (k^2 + 1) / (k^2 - 1)

# Now, we substitute the value of k=2 into the equation.

# Calculate the numerator of the formula
numerator = k**2 + 1

# Calculate the denominator of the formula
denominator = k**2 - 1

# Calculate the final ratio as a floating-point number
result = numerator / denominator

print("The formula for the ratio of maximum tangential stress to internal pressure is:")
print("Ratio = (k^2 + 1) / (k^2 - 1), where k = r_o / r_i\n")

print("Substituting the given value k = 2:")
print(f"Ratio = ( {k}^2 + 1 ) / ( {k}^2 - 1 )")
print(f"Ratio = ( {k*k} + 1 ) / ( {k*k} - 1 )")
print(f"Ratio = {int(numerator)} / {int(denominator)}")

print(f"\nThe final ratio is 5/3, which is approximately {result:.4f}.")