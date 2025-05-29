import math

# Given values
radius = 6  # inches
thickness = 1/3  # inches
number_of_pieces = 12

# Calculate the volume of the entire pizza
total_volume = math.pi * (radius ** 2) * thickness

# Calculate the volume of one piece
volume_of_one_piece = total_volume / number_of_pieces

print(volume_of_one_piece)