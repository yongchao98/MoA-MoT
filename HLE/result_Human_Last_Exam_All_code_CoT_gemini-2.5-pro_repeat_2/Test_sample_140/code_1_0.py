import math

def calculate_blocks_in_region(region_dims, block_dims):
    """Calculates the maximum number of blocks that fit in a region by checking all 3 orientations."""
    rx, ry, rz = sorted(region_dims, reverse=True)
    bx, by, bz = sorted(block_dims, reverse=True) # (9, 1, 1)

    # Only need to check one orientation since we sorted the dimensions.
    # The longest side of the block must align with a side of the region that is long enough.
    if rx < bx:
        return 0
    
    # Pack along the longest dimension
    fit = math.floor(rx / bx) * math.floor(ry / by) * math.floor(rz / bz)
    return fit

# Define dimensions
cube_dim = 15
block_dims = (9, 1, 1)
split_dim = block_dims[0]
remaining_dim = cube_dim - split_dim

# Define the 8 sub-regions based on the 9+6 split
sub_regions = {
    "9x9x9": (split_dim, split_dim, split_dim),
    "9x9x6": (split_dim, split_dim, remaining_dim),
    "9x6x9": (split_dim, remaining_dim, split_dim),
    "6x9x9": (remaining_dim, split_dim, split_dim),
    "9x6x6": (split_dim, remaining_dim, remaining_dim),
    "6x9x6": (remaining_dim, split_dim, remaining_dim),
    "6x6x9": (remaining_dim, remaining_dim, split_dim),
    "6x6x6": (remaining_dim, remaining_dim, remaining_dim)
}

total_blocks = 0
calculation_parts = []

print("Calculating the number of blocks for each sub-region:")
for name, dims in sub_regions.items():
    num_blocks = calculate_blocks_in_region(dims, block_dims)
    print(f"- Region {name}: {num_blocks} blocks")
    if num_blocks > 0:
        calculation_parts.append(str(num_blocks))
    total_blocks += num_blocks

# Final equation output
final_equation = " + ".join(calculation_parts)
print(f"\nTotal blocks = {final_equation} = {total_blocks}")
<<<351>>>