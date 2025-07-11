import collections

def solve_grid_assignments():
    """
    Calculates the number of valid 0/1 assignments for a 4x4 grid.
    """
    n = 4  # Number of rows
    m = 4  # Number of columns

    # Step 1: Find all valid row configurations (masks) for a row of length m
    # A mask is valid if no two adjacent bits are 1.
    valid_masks = []
    for i in range(1 << m):
        if (i & (i << 1)) == 0:
            valid_masks.append(i)
    
    # Sort masks for consistent output order
    valid_masks.sort()

    # Step 2: Base case - DP for the first row (n=1)
    # For the first row, each valid mask configuration is possible in one way.
    dp = {mask: 1 for mask in valid_masks}
    
    print(f"For a 1x4 grid, the total number of valid assignments is {sum(dp.values())}.")
    
    # Pre-compute compatibilities to speed up the main loop
    # Two masks are compatible if their bitwise AND is 0.
    compat_map = collections.defaultdict(list)
    for m1 in valid_masks:
        for m2 in valid_masks:
            if (m1 & m2) == 0:
                compat_map[m1].append(m2)

    # Step 3: Iterate from the second row to the n-th row
    for i in range(1, n):
        current_row_num = i + 1
        next_dp = {mask: 0 for mask in valid_masks}
        
        for mask in valid_masks:
            # Sum up the counts from compatible previous row masks
            count = sum(dp[prev_mask] for prev_mask in compat_map[mask])
            next_dp[mask] = count
            
        dp = next_dp
        total_assignments = sum(dp.values())
        print(f"For a {current_row_num}x4 grid, the total number of valid assignments is {total_assignments}.")

    # Step 4: Final calculation for the 4x4 grid
    final_counts = [dp[mask] for mask in valid_masks]
    final_total = sum(final_counts)

    print("\nTo find the total for a 4x4 grid, we sum the possibilities for each valid configuration of the last row:")
    
    # Display the final equation
    equation_str = " + ".join(map(str, final_counts))
    print(f"Total = {equation_str}")
    print(f"Total = {final_total}")


solve_grid_assignments()
<<<1157>>>