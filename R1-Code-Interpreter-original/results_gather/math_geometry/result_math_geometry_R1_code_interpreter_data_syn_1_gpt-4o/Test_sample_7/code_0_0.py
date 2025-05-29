import math

# Given values
diameter = 12
thickness = 1/3
number_of_pieces = 12

# Calculate radius
radius = diameter / 2

# Calculate total volume of the pizza
total_volume = math.pi * (radius ** 2) * thickness

# Calculate volume of one piece
volume_of_one_piece = total_volume / number_of_pieces

# Output the result
print(volume_of_one_piece)