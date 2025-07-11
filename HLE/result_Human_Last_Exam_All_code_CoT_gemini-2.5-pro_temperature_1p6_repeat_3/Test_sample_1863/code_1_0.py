import math

# Define the diameter of the full sphere from which the package is made.
package_diameter = 250.0

# Using the derived geometric relationship: d = D / (sqrt(3) + 1)
# where 'd' is the inner sphere's diameter and 'D' is the package's original sphere diameter.
sqrt_3 = math.sqrt(3)

# Calculate the theoretical maximum diameter of the sphere that can fit.
theoretical_max_diameter = package_diameter / (sqrt_3 + 1)

# The available spheres have diameters in 0.01 cm increments.
# We find the largest multiple of 0.01 that is less than or equal to the theoretical maximum.
# This is achieved by taking the floor of the value after scaling by 100.
final_diameter = math.floor(theoretical_max_diameter * 100) / 100

print("This script calculates the maximum diameter of a sphere that fits in a quarter-sphere package.")
print("\n--- Calculation Breakdown ---")
print(f"The calculation is based on the equation: d = D / (sqrt(3) + 1)")
print(f"Package's original sphere diameter (D): {package_diameter}")
print(f"Value of sqrt(3): {sqrt_3}")
print(f"The full equation with substituted numbers is:")
print(f"d = {package_diameter} / ({sqrt_3} + 1)")
print(f"Theoretical maximum diameter (d) \u2248 {theoretical_max_diameter} cm")

print("\n--- Final Answer ---")
print("Considering the 0.01 cm increments, we must truncate the theoretical diameter.")
print(f"Maximum allowed diameter = floor({theoretical_max_diameter} * 100) / 100")
print(f"The maximum diameter of a sphere that can fit in the package is {final_diameter:.2f} cm.")
