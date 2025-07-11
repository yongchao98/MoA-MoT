import math

def calculate_max_blocks(prism_dims, block_dims):
    """Calculates the maximum number of blocks that can fit in a prism by checking all orientations."""
    l, w, h = prism_dims
    bl, bw, bh = block_dims # Block length, width, height

    orientations = [
        (bl, bw, bh),
        (bl, bh, bw),
        (bw, bl, bh),
        (bw, bh, bl),
        (bh, bl, bw),
        (bh, bw, bl),
    ]

    max_blocks = 0
    # Check each possible orientation of the block
    for obl, obw, obh in set(orientations):
        # Check 3 possible alignments of the oriented block within the prism
        # Alignment 1
        if l >= obl and w >= obw and h >= obh:
             max_blocks = max(max_blocks, math.floor(l/obl) * math.floor(w/obw) * math.floor(h/obh))
        # Alignment 2 (swapping width and height)
        if l >= obl and w >= obh and h >= obw:
             max_blocks = max(max_blocks, math.floor(l/obl) * math.floor(w/obh) * math.floor(h/obw))
        # etc... but since two block dimensions are 1, we only need to check the 3 main axes
    
    # Simplified check for a Lx1x1 block
    o1 = math.floor(l / bl) * math.floor(w / bw) * math.floor(h / bh)
    o2 = math.floor(l / bw) * math.floor(w / bl) * math.floor(h / bh)
    o3 = math.floor(l / bw) * math.floor(w / bh) * math.floor(h / bl)
    
    return max(o1, o2, o3)


# Block and Cube dimensions
block_dims = (9, 1, 1)

# Decompose the cube and calculate blocks for each part
# Part 1: 9x9x9 cube
blocks_999 = calculate_max_blocks((9, 9, 9), block_dims)

# Part 2: 9x9x6 prism
blocks_996 = calculate_max_blocks((9, 9, 6), block_dims)

# Part 3: 9x6x6 prism
blocks_966 = calculate_max_blocks((9, 6, 6), block_dims)

# Part 4: 6x6x6 cube
blocks_666 = calculate_max_blocks((6, 6, 6), block_dims)

# There is 1 cube of 9x9x9, 3 prisms of 9x9x6, 3 prisms of 9x6x6, and 1 cube of 6x6x6
total_blocks = (1 * blocks_999) + (3 * blocks_996) + (3 * blocks_966) + (1 * blocks_666)

print("Calculating the total number of blocks based on cube decomposition:")
print(f"1 x (blocks in 9x9x9) = 1 * {blocks_999}")
print(f"3 x (blocks in 9x9x6) = 3 * {blocks_996}")
print(f"3 x (blocks in 9x6x6) = 3 * {blocks_966}")
print(f"1 x (blocks in 6x6x6) = 1 * {blocks_666}")
print("-" * 20)
# Final Equation
print(f"Total = (1 * {blocks_999}) + (3 * {blocks_996}) + (3 * {blocks_966}) + (1 * {blocks_666}) = {total_blocks}")
print(f"\nThe largest number of blocks that will fit is {total_blocks}.")
