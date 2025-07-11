def count_blocks(sequence: str) -> tuple[int, int]:
    """
    Counts the number of 0-blocks and 1-blocks in a binary string.
    
    A block is a sequence of one or more consecutive identical digits.
    For example, "00100011" has 2 blocks of 0s and 2 blocks of 1s.
    
    Returns:
        A tuple (num_0_blocks, num_1_blocks).
    """
    if not sequence:
        return 0, 0
    
    num_0_blocks = 0
    num_1_blocks = 0
    
    # Count the first block
    if sequence[0] == '0':
        num_0_blocks = 1
    else:
        num_1_blocks = 1
        
    # Iterate through the rest of the sequence to find new blocks
    for i in range(1, len(sequence)):
        # A new block starts if the current digit is different from the previous one
        if sequence[i] != sequence[i-1]:
            if sequence[i] == '0':
                num_0_blocks += 1
            else:
                num_1_blocks += 1
                
    return num_0_blocks, num_1_blocks

def solve_worst_case():
    """
    Calculates the minimum number of operations n needed for any transformation
    by finding the cost for the worst-case scenario.
    """
    # The worst-case sequences are those with the maximum number of blocks.
    # A 100-digit alternating sequence like "0101..." or "1010..." has
    # the maximum number of blocks: 50 of 0s and 50 of 1s.
    s1_worst = "01" * 50
    s2_worst = "10" * 50
    
    s1_0_blocks, s1_1_blocks = count_blocks(s1_worst)
    s2_0_blocks, s2_1_blocks = count_blocks(s2_worst)
    
    # Cost to transform via an all-zeros sequence
    # (remove all 1-blocks from s1, then create all 1-blocks for s2)
    cost_path_via_zeros = s1_1_blocks + s2_1_blocks
    
    # Cost to transform via an all-ones sequence
    # (remove all 0-blocks from s1, then create all 0-blocks for s2)
    cost_path_via_ones = s1_0_blocks + s2_0_blocks
    
    # The minimum operations required is the minimum of the two path costs
    n = min(cost_path_via_zeros, cost_path_via_ones)
    
    print("To find n, we calculate the cost for the worst-case sequences:")
    print(f"Worst-case S1: '{s1_worst[:10]}...'")
    print(f"Worst-case S2: '{s2_worst[:10]}...'")
    print("-" * 20)
    print(f"Number of 0-blocks in S1: {s1_0_blocks}")
    print(f"Number of 1-blocks in S1: {s1_1_blocks}")
    print(f"Number of 0-blocks in S2: {s2_0_blocks}")
    print(f"Number of 1-blocks in S2: {s2_1_blocks}")
    print("-" * 20)
    
    print("Path 1 cost (transforming via all '0's):")
    print(f"num_1_blocks(S1) + num_1_blocks(S2) = {s1_1_blocks} + {s2_1_blocks} = {cost_path_via_zeros}")
    
    print("\nPath 2 cost (transforming via all '1's):")
    print(f"num_0_blocks(S1) + num_0_blocks(S2) = {s1_0_blocks} + {s2_0_blocks} = {cost_path_via_ones}")
    
    print("\nThe maximum number of operations 'n' is the minimum of these two paths for the worst case:")
    print(f"n = min({cost_path_via_zeros}, {cost_path_via_ones}) = {n}")

solve_worst_case()
