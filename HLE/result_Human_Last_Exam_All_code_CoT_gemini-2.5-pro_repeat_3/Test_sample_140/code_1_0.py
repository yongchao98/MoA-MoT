import math

# Define the dimensions
cube_edge = 15
block_length = 9
block_width = 1
block_height = 1

# Step 1: Fill a 9x15x15 section of the cube.
# The 9-unit length of the blocks is aligned with the first 15-unit edge of the cube.
num1 = (block_length // block_length) * (cube_edge // block_width) * (cube_edge // block_height)

# The remaining space has dimensions (15-9)x15x15 = 6x15x15
remaining_dim_1 = cube_edge - block_length

# Step 2: Fill the 6x15x15 space.
# We align the 9-unit length of the blocks with the second 15-unit edge.
# The number of blocks is calculated for a 6x15x15 space, with blocks oriented as 1x9x1.
num2 = (remaining_dim_1 // block_width) * (cube_edge // block_length) * (cube_edge // block_height)

# The remaining space has dimensions 6x(15-9)x15 = 6x6x15
remaining_dim_2 = cube_edge - block_length

# Step 3: Fill the 6x6x15 space.
# We align the 9-unit length of the blocks with the last 15-unit edge.
# The number of blocks is calculated for a 6x6x15 space, with blocks oriented as 1x1x9.
num3 = (remaining_dim_1 // block_width) * (remaining_dim_2 // block_height) * (cube_edge // block_length)

# The final remaining space is 6x6x(15-9) = 6x6x6, where no more blocks can fit.

# Step 4: Calculate the total and print the final equation.
total_blocks = num1 + num2 + num3

print("The problem can be solved by splitting the cube into sections:")
print(f"Section 1 (9x15x15): We can fit {num1} blocks.")
print(f"Section 2 (6x15x15): We can fit {num2} blocks by changing orientation.")
print(f"Section 3 (6x6x15): We can fit {num3} blocks with a third orientation.")
print("The final remaining 6x6x6 space cannot fit any 9x1x1 blocks.")
print("\nThe largest number of blocks is the sum of the blocks from all sections:")
print(f"{num1} + {num2} + {num3} = {total_blocks}")
