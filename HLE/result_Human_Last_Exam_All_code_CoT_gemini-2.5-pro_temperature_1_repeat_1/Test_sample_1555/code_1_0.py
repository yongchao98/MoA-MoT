import math

# This script calculates the leading-order asymptotic mass ratio between the
# lightest and the next-to-lightest physical excitations in the CP(N-1) model.

# The exact mass ratio between the k=2 and k=1 states is:
# Ratio(N) = m_2 / m_1 = 2 * cos(pi / N)

# We need to find the asymptotic value of this ratio as N approaches infinity.
# This involves calculating the limit of the expression as N -> infinity.
# In this limit, the term (pi / N) approaches 0.

# Define the terms for the final calculation based on the limit.
# The constant factor in the expression.
factor = 2
# The argument of the cosine function in the limit (pi/N -> 0).
limit_arg = 0

# Calculate the value of cos(0).
cos_value = math.cos(limit_arg)

# Calculate the final result of the limit.
asymptotic_ratio = factor * cos_value

print("The asymptotic mass ratio is the limit of 2 * cos(pi / N) as N -> infinity.")
print("The final calculation is:")
# The f-string will print the equation with the evaluated numbers.
print(f"{factor} * cos({limit_arg}) = {asymptotic_ratio}")