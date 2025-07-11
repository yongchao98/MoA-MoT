import math

def solve_packing():
    """
    Calculates the largest number of 9x1x1 blocks that fit in a 15x15x15 cube.
    """

    # --- Problem Dimensions ---
    cube_edge = 15
    block_len, block_w, block_h = 9, 1, 1

    # Stage 1: Pack the initial 15x15x15 cube
    # Orient blocks as 9x1x1 along one axis.
    # We fit a layer of blocks that is floor(15/9) = 1 block deep.
    # The layer is 15 blocks wide and 15 blocks high.
    num_blocks_1 = math.floor(15 / block_len) * math.floor(15 / block_w) * math.floor(15 / block_h)
    
    # This leaves a 6x15x15 space (15-9=6).
    
    # Stage 2: Pack the remaining 6x15x15 space
    # We must reorient the blocks. Place the 9-unit side along a 15-unit axis.
    # We can fit floor(15/9)=1 column of blocks.
    # The cross-section of this column is 6x15 blocks.
    num_blocks_2 = math.floor(6 / block_w) * math.floor(15 / block_len) * math.floor(15 / block_h)

    # This leaves a 6x6x15 space (15-9=6).

    # Stage 3: Pack the remaining 6x6x15 space
    # We again reorient the blocks. Place the 9-unit side along the last 15-unit axis.
    # We can fit floor(15/9)=1 set of blocks deep.
    # The face of this set is 6x6 blocks.
    num_blocks_3 = math.floor(6 / block_w) * math.floor(6 / block_h) * math.floor(15 / block_len)
    
    # This leaves a 6x6x6 space, into which no 9-unit block can fit.
    
    # Calculate the total
    total_blocks = num_blocks_1 + num_blocks_2 + num_blocks_3
    
    print("The total number of blocks is the sum of blocks from each stage of packing.")
    print(f"Final Calculation: {num_blocks_1} + {num_blocks_2} + {num_blocks_3} = {total_blocks}")

solve_packing()