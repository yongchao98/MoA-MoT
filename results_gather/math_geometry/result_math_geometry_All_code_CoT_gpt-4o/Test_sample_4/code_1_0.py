# Given values
cube_edge_length = 4
pyramid_base_edge_length = 8

# Calculate the volume of the cube
cube_volume = cube_edge_length ** 3

# Calculate the volume of the pyramid in terms of h
# Volume of pyramid = (1/3) * base_area * height
# Base area = pyramid_base_edge_length ** 2
pyramid_volume = (1/3) * (pyramid_base_edge_length ** 2) * 1  # Volume per unit height

# Set the volumes equal and solve for h
h = cube_volume / pyramid_volume

# Output the result
print(f"<<<{h}>>>")