# Calculate the volume of the sphere
radius_sphere = 5
volume_sphere = (4/3) * 3.141592653589793 * (radius_sphere ** 3)

# Calculate the volume of the cylinder
radius_cylinder = 3
height_cylinder = 8
volume_cylinder = 3.141592653589793 * (radius_cylinder ** 2) * height_cylinder

# Calculate the volume of the space outside the cylinder but inside the sphere
volume_space = volume_sphere - volume_cylinder

# Calculate W
W = volume_space / 3.141592653589793

print(W)