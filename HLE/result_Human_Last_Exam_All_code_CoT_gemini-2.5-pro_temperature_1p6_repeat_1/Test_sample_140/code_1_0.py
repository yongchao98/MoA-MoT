import math

def calculate_packing():
    """
    Calculates the maximum number of 9x1x1 blocks that fit in a 15x15x15 cube.
    """
    print("This problem is solved by partitioning the 15x15x15 cube for optimal packing.")
    print("A side of length 15 is best split into a 9-unit part and a 6-unit part.")
    print("Applying this to all three axes divides the cube into 8 smaller cuboids.")
    print("We will now calculate the blocks for each part and sum them up.\n")

    # --- Calculation for the 9x9x9 cuboid ---
    # There is one cuboid of this size.
    # Blocks can be oriented along any axis.
    # Using integer division to find how many blocks fit.
    blocks_9x9x9 = (9 // 9) * (9 // 1) * (9 // 1)
    print("1. For the single 9x9x9 cuboid:")
    print(f"   Calculation: ({9}//{9}) * ({9}//{1}) * ({9}//{1}) = {blocks_9x9x9} blocks")
    print("-" * 20)

    # --- Calculation for the three 9x9x6 cuboids ---
    # To fit a 9-unit block, it must be aligned with a 9-unit side.
    blocks_per_9x9x6 = (9 // 9) * (9 // 1) * (6 // 1)
    num_9x9x6 = 3
    total_blocks_9x9x6 = num_9x9x6 * blocks_per_9x9x6
    print(f"2. For the {num_9x9x6} cuboids of size 9x9x6:")
    print(f"   Blocks per cuboid: ({9}//{9}) * ({9}//{1}) * ({6}//{1}) = {blocks_per_9x9x6} blocks")
    print(f"   Subtotal: {num_9x9x6} * {blocks_per_9x9x6} = {total_blocks_9x9x6} blocks")
    print("-" * 20)

    # --- Calculation for the three 9x6x6 cuboids ---
    # The 9-unit block must be aligned with the single 9-unit side.
    blocks_per_9x6x6 = (9 // 9) * (6 // 1) * (6 // 1)
    num_9x6x6 = 3
    total_blocks_9x6x6 = num_9x6x6 * blocks_per_9x6x6
    print(f"3. For the {num_9x6x6} cuboids of size 9x6x6:")
    print(f"   Blocks per cuboid: ({9}//{9}) * ({6}//{1}) * ({6}//{1}) = {blocks_per_9x6x6} blocks")
    print(f"   Subtotal: {num_9x6x6} * {blocks_per_9x6x6} = {total_blocks_9x6x6} blocks")
    print("-" * 20)
    
    # --- Calculation for the 6x6x6 cuboid ---
    # The longest side is 6, so a 9-unit block cannot fit.
    blocks_6x6x6 = 0
    print("4. For the single 6x6x6 cuboid:")
    print(f"   The longest side is 6, which is less than the block length of 9.")
    print(f"   Calculation: {blocks_6x6x6} blocks")
    print("-" * 20)

    # --- Final Summation ---
    total_blocks = blocks_9x9x9 + total_blocks_9x9x6 + total_blocks_9x6x6 + blocks_6x6x6
    print("\nThe total number of blocks is the sum from all parts.")
    print(f"Final Equation: {blocks_9x9x9} + {total_blocks_9x9x6} + {total_blocks_9x6x6} + {blocks_6x6x6}")
    print(f"Largest number of blocks = {total_blocks}")

calculate_packing()