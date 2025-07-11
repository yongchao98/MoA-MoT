import math

# Define the physical constants and given values
# Bohr radius (radius of the first orbit), r1, in meters
r1 = 5.29177e-11  # meters
# Uncertainty in position, delta_x, in meters
delta_x_pm = 10  # picometers
delta_x_m = delta_x_pm * 1e-12  # convert pm to m

# The ratio of momentum uncertainty to momentum is given by r1 / (2 * delta_x)
# Let's calculate the components of the final equation
numerator = r1
denominator = 2 * delta_x_m

# Calculate the final ratio
ratio = numerator / denominator

# Print the final equation with all the numbers as requested
print("The ratio of the uncertainty of the momentum to the momentum is given by the formula:")
print("Ratio = r₁ / (2 * Δx)")
print(f"Ratio = {numerator:.3e} m / (2 * {delta_x_m:.3e} m)")
print(f"Ratio = {numerator:.3e} m / {denominator:.3e} m")
print(f"Final Ratio = {ratio}")
