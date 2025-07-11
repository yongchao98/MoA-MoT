import math

def solve_packing_problem():
    """
    Calculates the maximum number of 9x1x1 blocks in a 15x15x15 cube.
    """
    block_l, block_w, block_h = 9, 1, 1

    # This function calculates the max blocks for a given region by checking 3 orientations.
    def calculate_max_blocks_in_region(l, w, h):
        # Orientation 1: 9x1x1
        fit1 = math.floor(l / 9) * math.floor(w / 1) * math.floor(h / 1)
        # Orientation 2: 1x9x1
        fit2 = math.floor(l / 1) * math.floor(w / 9) * math.floor(h / 1)
        # Orientation 3: 1x1x9
        fit3 = math.floor(l / 1) * math.floor(w / 1) * math.floor(h / 9)
        return max(fit1, fit2, fit3)

    print("The 15x15x15 cube is partitioned into 8 smaller regions to pack the 9x1x1 blocks efficiently.")
    print("We calculate the number of blocks that fit in each type of region:\n")

    # 1. The 9x9x9 core region
    blocks_999 = calculate_max_blocks_in_region(9, 9, 9)
    print(f"In the single 9x9x9 region, we can fit:")
    print(f"  floor(9/9) * floor(9/1) * floor(9/1) = 1 * 9 * 9 = {blocks_999} blocks.\n")

    # 2. The three 9x9x6 regions
    blocks_996 = calculate_max_blocks_in_region(9, 9, 6)
    total_996 = 3 * blocks_996
    print(f"In the three 9x9x6 regions, we can fit:")
    print(f"  floor(9/9) * floor(9/1) * floor(6/1) = 1 * 9 * 6 = {blocks_996} blocks per region.")
    print(f"  Subtotal for these 3 regions = 3 * {blocks_996} = {total_996} blocks.\n")

    # 3. The three 9x6x6 regions
    blocks_966 = calculate_max_blocks_in_region(9, 6, 6)
    total_966 = 3 * blocks_966
    print(f"In the three 9x6x6 regions, we can fit:")
    print(f"  floor(9/9) * floor(6/1) * floor(6/1) = 1 * 6 * 6 = {blocks_966} blocks per region.")
    print(f"  Subtotal for these 3 regions = 3 * {blocks_966} = {total_966} blocks.\n")

    # 4. The 6x6x6 corner region
    blocks_666 = calculate_max_blocks_in_region(6, 6, 6)
    print(f"In the single 6x6x6 region, all dimensions are smaller than 9, so we can fit:")
    print(f"  {blocks_666} blocks.\n")

    # 5. Final Sum
    total_blocks = blocks_999 + total_996 + total_966 + blocks_666
    
    print("The largest number of blocks is the sum from all regions:")
    # The final equation shows the contribution from each group of regions
    print(f"Total = {blocks_999} + {total_996} + {total_966} + {blocks_666} = {total_blocks}")

solve_packing_problem()