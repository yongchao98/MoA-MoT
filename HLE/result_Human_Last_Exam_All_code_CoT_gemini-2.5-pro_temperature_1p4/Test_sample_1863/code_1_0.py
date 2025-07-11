import math

# The diameter of the original sphere from which the quarter-sphere is made.
package_diameter = 250.0

# The constant values in our derived formula.
val_sqrt = 3
val_add = 1

# The formula for the diameter of the inscribed sphere is:
# Diameter = package_diameter / (sqrt(3) + 1)
max_sphere_diameter = package_diameter / (math.sqrt(val_sqrt) + val_add)

# Output the equation with all the numbers and the final result.
print("To find the maximum diameter of the sphere that fits in the package, we use the following equation:")
print(f"{package_diameter} / (sqrt({val_sqrt}) + {val_add})")
print("\nResult:")
print(f"The maximum diameter of the sphere is {max_sphere_diameter} cm.")
