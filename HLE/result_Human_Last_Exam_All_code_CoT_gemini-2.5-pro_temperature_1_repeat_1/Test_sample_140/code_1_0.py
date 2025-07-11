import math

def solve_packing_problem():
    """
    Calculates the maximum number of 9x1x1 blocks that can fit in a 15x15x15 cube.
    """
    cube_edge = 15
    block_l, block_w, block_h = 9, 1, 1

    print(f"This program calculates the largest number of {block_l}x{block_w}x{block_h} blocks that fit inside a {cube_edge}x{cube_edge}x{cube_edge} cube.")
    print("The strategy is to partition the cube based on the block's longest dimension.")
    print(f"We partition each 15-unit edge into a {block_l}-unit part and a {cube_edge - block_l}-unit part (15 = 9 + 6).\n")

    # Partition 1: 9x9x9 cuboid
    # There is 1 of these.
    p1_dims = (9, 9, 9)
    num_p1 = 1
    # Blocks must be oriented as 9x1x1.
    blocks_p1 = (p1_dims[0] // block_l) * (p1_dims[1] // block_w) * (p1_dims[2] // block_h)
    total_blocks_p1 = num_p1 * blocks_p1
    print(f"1. In the single {p1_dims[0]}x{p1_dims[1]}x{p1_dims[2]} cuboid, we can fit:")
    print(f"   ({p1_dims[0]} // {block_l}) * ({p1_dims[1]} // {block_w}) * ({p1_dims[2]} // {block_h}) = {blocks_p1} blocks.")

    # Partition 2: 9x9x6 cuboids
    # There are 3 of these (9x9x6, 9x6x9, 6x9x9).
    p2_dims = (9, 9, 6)
    num_p2 = 3
    # Blocks must be oriented along a 9-unit edge.
    blocks_p2 = (p2_dims[0] // block_l) * (p2_dims[1] // block_w) * (p2_dims[2] // block_h)
    total_blocks_p2 = num_p2 * blocks_p2
    print(f"\n2. In the three {p2_dims[0]}x{p2_dims[1]}x{p2_dims[2]} cuboids, we can fit:")
    print(f"   {num_p2} * (({p2_dims[0]} // {block_l}) * ({p2_dims[1]} // {block_w}) * ({p2_dims[2]} // {block_h})) = {num_p2} * {blocks_p2} = {total_blocks_p2} blocks.")

    # Partition 3: 9x6x6 cuboids
    # There are 3 of these (9x6x6, 6x9x6, 6x6x9).
    p3_dims = (9, 6, 6)
    num_p3 = 3
    # Blocks must be oriented along the 9-unit edge.
    blocks_p3 = (p3_dims[0] // block_l) * (p3_dims[1] // block_w) * (p3_dims[2] // block_h)
    total_blocks_p3 = num_p3 * blocks_p3
    print(f"\n3. In the three {p3_dims[0]}x{p3_dims[1]}x{p3_dims[2]} cuboids, we can fit:")
    print(f"   {num_p3} * (({p3_dims[0]} // {block_l}) * ({p3_dims[1]} // {block_w}) * ({p3_dims[2]} // {block_h})) = {num_p3} * {blocks_p3} = {total_blocks_p3} blocks.")

    # Partition 4: 6x6x6 cuboid
    # There is 1 of these.
    p4_dims = (6, 6, 6)
    total_blocks_p4 = 0
    print(f"\n4. In the single {p4_dims[0]}x{p4_dims[1]}x{p4_dims[2]} cuboid, we can fit:")
    print(f"   {total_blocks_p4} blocks, because all its dimensions are smaller than {block_l}.")

    # Final calculation
    total_blocks = total_blocks_p1 + total_blocks_p2 + total_blocks_p3 + total_blocks_p4
    print("\n-----------------------------------------------------")
    print("The total number of blocks is the sum from all partitions:")
    print(f"Total = {total_blocks_p1} + {total_blocks_p2} + {total_blocks_p3} + {total_blocks_p4} = {total_blocks}")
    print("-----------------------------------------------------")

solve_packing_problem()
<<<351>>>