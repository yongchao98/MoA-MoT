import math

# The diameter of the quarter-sphere package is given.
package_diameter = 250.0

# The problem is to find the maximum diameter of a sphere that can be inscribed
# in a quarter-sphere. Let 'D' be the diameter of the quarter-sphere and 'd' be
# the diameter of the inscribed sphere.
# The geometric relationship for the largest possible inscribed sphere is:
# d = D / (sqrt(2) + 1)

# We will now calculate this value.
# First, calculate the components of the equation.
val_sqrt_2 = math.sqrt(2)
denominator = val_sqrt_2 + 1

# Now, calculate the final diameter 'd'.
max_sphere_diameter = package_diameter / denominator

# Display the final equation with the numbers substituted.
print("To find the maximum diameter 'd' of a sphere that fits in a quarter-sphere of diameter D = {} cm, we use the formula:".format(package_diameter))
print("\nd = D / (sqrt(2) + 1)")
print("\nSubstituting the numbers into the equation:")
print("d = {} / ({} + 1)".format(package_diameter, val_sqrt_2))
print("d = {} / {}".format(package_diameter, denominator))
print("\nCalculation Result:")
print("The maximum diameter of the sphere is {} cm.".format(max_sphere_diameter))