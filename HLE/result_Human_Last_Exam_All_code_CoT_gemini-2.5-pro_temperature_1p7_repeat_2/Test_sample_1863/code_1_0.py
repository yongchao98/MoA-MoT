import math

# The diameter of the quarter-sphere package is given.
package_diameter = 250.0

# We use the derived formula to find the maximum diameter (d) of the sphere that can fit inside.
# The formula is: d = D / (sqrt(3) + 1), where D is the package diameter.

# Calculate the value of the denominator.
denominator = math.sqrt(3) + 1

# Calculate the maximum diameter of the inscribed sphere.
sphere_diameter = package_diameter / denominator

# Print the final equation with each number and the result.
print("The formula for the maximum diameter (d) of a sphere inside a quarter-sphere of diameter (D) is:")
print("d = D / (sqrt(3) + 1)")
print("\nSubstituting the given values:")
print(f"d = {package_diameter} / (√3 + 1)")
print(f"d = {package_diameter} / ({denominator})")
print(f"\nAfter performing the calculation, the final answer is:")
print(f"The maximum diameter of the sphere is {sphere_diameter} cm.")

print(f"\n<<<{sphere_diameter}>>>")