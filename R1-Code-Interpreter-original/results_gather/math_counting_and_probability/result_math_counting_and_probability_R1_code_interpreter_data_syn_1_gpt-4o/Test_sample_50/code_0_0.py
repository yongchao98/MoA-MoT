# Initial number of tiles
initial_blue_tiles = 13
initial_green_tiles = 6

# Side length of the original hexagon
original_side_length = 2

# Perimeter of the original hexagon
original_perimeter = 6 * original_side_length

# New side length after adding a border
new_side_length = original_side_length + 1

# Perimeter of the new hexagon
new_perimeter = 6 * new_side_length

# Additional green tiles added by the border
additional_green_tiles = new_perimeter - original_perimeter

# Total green tiles after adding the border
total_green_tiles = initial_green_tiles + additional_green_tiles

# Difference between total green tiles and total blue tiles
difference = total_green_tiles - initial_blue_tiles

print(difference)