# Billet dimensions in cm
billet_length = 16
billet_width = 11
billet_height = 4

# T1 cube dimensions and price
t1_side = 1
t1_price = 5

# Calculate the maximum number of T1 cubes that can fit in the billet
# This is equivalent to dividing the billet volume by the cube volume
num_t1_cubes_length = billet_length // t1_side
num_t1_cubes_width = billet_width // t1_side
num_t1_cubes_height = billet_height // t1_side

total_t1_cubes = num_t1_cubes_length * num_t1_cubes_width * num_t1_cubes_height

# Calculate the maximum total value
max_value = total_t1_cubes * t1_price

# Print the final equation as requested
print(f"The optimal solution is to fill the billet with T1 cubes.")
print(f"Number of cubes: {num_t1_cubes_length} * {num_t1_cubes_width} * {num_t1_cubes_height} = {total_t1_cubes}")
print(f"Total value equation: {total_t1_cubes} * {t1_price} = {max_value}")
