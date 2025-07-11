import math

# Define the diameter of the quarter-sphere package
package_diameter = 250.0

# Calculate the radius of the package
package_radius = package_diameter / 2

# Calculate sqrt(3) for the formula
sqrt_3 = math.sqrt(3)

# The maximum theoretical diameter 'd' of the sphere that can fit is given by the formula:
# d = package_diameter / (sqrt(3) + 1)
max_theoretical_diameter = package_diameter / (sqrt_3 + 1)

# The available spheres have diameters in increments of 0.01 cm.
# We must find the largest diameter that is a multiple of 0.01 and does not exceed the theoretical maximum.
# This is equivalent to truncating the theoretical diameter to two decimal places.
final_diameter = math.floor(max_theoretical_diameter * 100) / 100

# Print the final equation with the numbers plugged in and the result.
print("The equation to find the maximum theoretical diameter (d) is:")
print(f"d = {package_diameter} / (sqrt(3) + 1)")
print("\nPlugging in the value for sqrt(3):")
print(f"d = {package_diameter} / ({sqrt_3} + 1)")
print(f"Theoretical maximum diameter ≈ {max_theoretical_diameter} cm")
print("\nRounding down to the nearest available size (in increments of 0.01 cm):")
print(f"Final Answer: {final_diameter:.2f} cm")