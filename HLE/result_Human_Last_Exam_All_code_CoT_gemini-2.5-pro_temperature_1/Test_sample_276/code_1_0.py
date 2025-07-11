import math

# This problem boils down to calculating the area of four specific regions
# in the corners of the unit square [0,1]x[0,1].
# The total area (which equals the probability) can be expressed analytically as:
# P = 24 - 72 * arcsin(1/6) - 2 * sqrt(35)
# The code below calculates this value.

# Calculate the components of the formula
val_asin = math.asin(1/6)
val_sqrt = math.sqrt(35)

# Calculate the probability
probability = 24 - 72 * val_asin - 2 * val_sqrt

# Output the result formatted to four significant digits
print(f"The equation for the probability P is:")
print(f"P = 24 - 72 * arcsin(1/6) - 2 * sqrt(35)")
print(f"Intermediate values:")
print(f"arcsin(1/6) = {val_asin}")
print(f"sqrt(35) = {val_sqrt}")
print(f"The calculated probability is approximately: {probability:.4g}")
