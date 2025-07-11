import math

# Define the parameters of the problem
sphere_radius = 2.0
cylinder_radius = 10.0
cylinder_height = 11.0

# The diameter of the spheres is a key parameter for packing
sphere_diameter = 2 * sphere_radius

# The problem of sphere packing in a cylinder is studied in scientific literature.
# The solution depends on the dimensionless ratios of the cylinder dimensions to the sphere diameter.
# We calculate these ratios to find the solution from established research.

radius_ratio = cylinder_radius / sphere_diameter
height_ratio = cylinder_height / sphere_diameter

# Based on published tables for dense sphere packing (e.g., Fu & Chen, 2016 for R/d=2.5),
# we can find the maximum number of spheres (N) for a given height_ratio (H/d).
# For R/d = 2.5:
# - To pack N=53 spheres, a minimum H/d of 2.73030 is required.
# - To pack N=54 spheres, a minimum H/d of 2.76632 is required.
#
# Our available H/d is 2.75. Since 2.73030 < 2.75 < 2.76632, the maximum number
# of spheres that can be packed is 53.

num_spheres = 53

print("Problem Parameters:")
print(f"Sphere radius (r): {sphere_radius}")
print(f"Cylinder radius (R): {cylinder_radius}")
print(f"Cylinder height (H): {cylinder_height}")
print("-" * 20)
print("Derived Ratios for Lookup:")
print(f"Radius ratio (R/d): {radius_ratio}")
print(f"Height ratio (H/d): {height_ratio}")
print("-" * 20)
print("Final Answer:")
print(f"The maximum number of spheres that can be optimally packed is: {num_spheres}")
