import math

# The diameter of the quarter-sphere package in cm.
package_diameter = 250.0

# The problem is to find the maximum diameter 'd' of a sphere that can fit inside
# a quarter-sphere of diameter 'D'.
# The relationship is given by the formula: d = D / (sqrt(3) + 1)

# Perform the calculation
val_sqrt_3 = math.sqrt(3)
denominator = val_sqrt_3 + 1
sphere_diameter = package_diameter / denominator

# Print the equation with the numbers substituted, as requested.
print(f"The equation to find the maximum sphere diameter (d) is:")
print(f"d = Package Diameter / (√3 + 1)")
print(f"Substituting the values:")
print(f"d = {package_diameter} / ({val_sqrt_3} + 1)")
print(f"d = {package_diameter} / {denominator}")
print(f"The final calculated maximum diameter is:")
print(f"d = {sphere_diameter}")