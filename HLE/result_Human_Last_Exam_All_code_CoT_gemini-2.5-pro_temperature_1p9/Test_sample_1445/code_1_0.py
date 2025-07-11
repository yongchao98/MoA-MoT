def solve_sequence_transformation():
    """
    Calculates the minimum number of operations n to transform any 100-digit binary
    sequence into any other.
    """
    
    # Let n(S, T) be the minimum operations to transform sequence S to T.
    # We will find an upper and a lower bound for n.

    # 1. Upper Bound Calculation
    # Any transformation S -> T can be done via an intermediate Z_0 (e.g., all zeros).
    # n(S, T) <= n(S, Z_0) + n(Z_0, T).
    # The cost n(S, Z_0) is the number of '1'-blocks in S, as we can remove them
    # one by one, causing the '0'-blocks to merge.
    
    # The maximum number of blocks of a single type in a 100-digit sequence
    # occurs in an alternating sequence like '1010...' or '0101...'.
    # For a length of 100, there can be at most 50 blocks of '1's.
    max_blocks_of_one_type = 100 // 2
    
    # So, n(S, Z_0) is at most 50. Since operations are reversible, n(Z_0, T) is also at most 50.
    cost_S_to_Z0 = max_blocks_of_one_type
    cost_Z0_to_T = max_blocks_of_one_type
    
    upper_bound = cost_S_to_Z0 + cost_Z0_to_T
    
    print("Upper bound calculation:")
    print(f"The maximum cost to transform S to Z_0 (all '0's) is {cost_S_to_Z0}.")
    print(f"The maximum cost to transform Z_0 to T is {cost_Z0_to_T}.")
    print(f"So, n <= n(S, Z_0) + n(Z_0, T) = {cost_S_to_Z0} + {cost_Z0_to_T} = {upper_bound}")

    # 2. Lower Bound Calculation
    # Consider a worst-case scenario: S = '0' * 100 and T = '01' * 50.
    # To transform S to T, we can:
    # a) Reduce S to a single '0'. This takes 1 operation (remove 99 '0's).
    # b) Build T from '0'. T has 100 blocks. We already have the first '0' block.
    #    We need to create the other 99 blocks by inserting '1', then '0', etc.
    #    This requires 99 insertion operations.
    cost_resize_S = 1
    cost_insertions_for_T = 99
    lower_bound = cost_resize_S + cost_insertions_for_T
    
    print("\nLower bound calculation:")
    print("Consider transforming S = '0'*100 to T = '01'*50.")
    print(f"Cost to reduce S to a single '0': {cost_resize_S}")
    print(f"Cost to build T from '0': {cost_insertions_for_T}")
    print(f"Total operations for this case = {cost_resize_S} + {cost_insertions_for_T} = {lower_bound}")

    # 3. Conclusion
    # Since the upper bound and lower bound are the same, this is our answer.
    print(f"\nThe upper bound ({upper_bound}) and lower bound ({lower_bound}) match.")
    print("\nFinal Answer:")
    print(f"The minimum number of operations n is {lower_bound}.")

solve_sequence_transformation()
<<<100>>>