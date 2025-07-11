# Set the spatial dimensionality 'd'.
# This value must be below the upper critical dimension of 4.
# We will use the common and physically relevant dimension d = 3.
d = 3

# In the epsilon expansion (where epsilon = 4 - d) for a scalar (N=1) field theory,
# the specific heat critical exponent, alpha, is given to the first order by the formula:
# alpha = (4 - d) / 6.

# Perform the calculation.
alpha_numerator = 4 - d
alpha_denominator = 6
alpha_value = alpha_numerator / alpha_denominator

# As requested, we will print the final equation showing each number involved in the calculation.
print(f"For a spatial dimensionality of d = {d}:")
print(f"The quantitative scaling relation is alpha = (4 - d) / 6.")
print(f"Plugging in the value for d, we get the equation:")
print(f"alpha = (4 - {d}) / 6 = {alpha_numerator} / {alpha_denominator} = {alpha_value:.4f}")