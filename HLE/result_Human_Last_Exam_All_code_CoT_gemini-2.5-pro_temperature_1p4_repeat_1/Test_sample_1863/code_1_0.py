import math

# Diameter of the quarter-sphere package in cm
package_diameter = 250

# The equation for the theoretical maximum diameter (d) of the inner sphere is:
# d = package_diameter / (sqrt(3) + 1)
theoretical_diameter = package_diameter / (math.sqrt(3) + 1)

# The available spheres have diameters in increments of 0.01 cm.
# We need to find the largest diameter that is less than or equal to the theoretical maximum.
# We can do this by flooring the value to two decimal places.
final_diameter = math.floor(theoretical_diameter * 100) / 100

# Print the final equation with the numbers and the result
print(f"The equation to find the maximum theoretical diameter (d) is:")
print(f"d = {package_diameter} / (sqrt(3) + 1)")
print(f"Calculated theoretical diameter: {theoretical_diameter:.4f} cm")
print(f"Given the available sphere increments, the maximum possible diameter is: {final_diameter:.2f} cm")
