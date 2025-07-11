# Define the dimensions of the cube and the block
cube_edge = 15
block_l, block_w, block_h = 9, 1, 1

# --- Step 1: First placement ---
# Fill a 9x15x15 section of the 15x15x15 cube.
# Blocks are oriented with their 9-unit length along one axis.
count1 = (cube_edge // block_l) * (cube_edge // block_w) * (cube_edge // block_h)
# This corresponds to 1 * 15 * 15 = 225 blocks.
print(f"In the first step, we place blocks in a {block_l}x{cube_edge}x{cube_edge} section.")
print(f"Number of blocks in this step: ({cube_edge}//{block_l}) * ({cube_edge}//{block_w}) * ({cube_edge}//{block_h}) = {count1}")
print("-" * 20)

# The remaining space is (15 - 9) x 15 x 15, which is 6x15x15.
rem_space1_x = cube_edge % block_l
rem_space1_y = cube_edge
rem_space1_z = cube_edge

# --- Step 2: Second placement in the remaining 6x15x15 space ---
# We orient the blocks differently, aligning the 9-unit side along the 15-unit edge.
count2 = (rem_space1_x // block_w) * (rem_space1_y // block_l) * (rem_space1_z // block_h)
# This corresponds to 6 * 1 * 15 = 90 blocks.
print(f"In the second step, we fill the remaining {rem_space1_x}x{rem_space1_y}x{rem_space1_z} space.")
print(f"Number of blocks in this step: ({rem_space1_x}//{block_w}) * ({rem_space1_y}//{block_l}) * ({rem_space1_z}//{block_h}) = {count2}")
print("-" * 20)

# The new remaining space is 6 x (15 - 9) x 15, which is 6x6x15.
rem_space2_x = rem_space1_x
rem_space2_y = rem_space1_y % block_l
rem_space2_z = rem_space1_z

# --- Step 3: Third placement in the remaining 6x6x15 space ---
# We orient the blocks again, aligning the 9-unit side along the last 15-unit edge.
count3 = (rem_space2_x // block_w) * (rem_space2_y // block_h) * (rem_space2_z // block_l)
# This corresponds to 6 * 6 * 1 = 36 blocks.
print(f"In the third step, we fill the remaining {rem_space2_x}x{rem_space2_y}x{rem_space2_z} space.")
print(f"Number of blocks in this step: ({rem_space2_x}//{block_w}) * ({rem_space2_y}//{block_h}) * ({rem_space2_z}//{block_l}) = {count3}")
print("-" * 20)

# The final remaining space is 6x6x6, which is too small for any more blocks.

# --- Total Calculation ---
total_blocks = count1 + count2 + count3
print("The final calculation is the sum of blocks from all steps:")
print(f"{count1} + {count2} + {count3} = {total_blocks}")
