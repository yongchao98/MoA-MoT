import collections

def solve_grid_assignments():
    """
    Calculates the number of valid 0/1 assignments for an n x m grid
    based on the given implication rules.
    """
    n = 4
    m = 4

    # Step 1: Find all valid row configurations (masks)
    # A mask is valid if no two adjacent bits are 1.
    valid_masks = []
    for mask in range(1 << m):
        if (mask & (mask << 1)) == 0:
            valid_masks.append(mask)

    # Step 2: Initialize DP table for the first row.
    # dp[mask] stores the number of ways to fill the grid up to the current row
    # with the current row having the configuration 'mask'.
    dp = collections.defaultdict(int)
    for mask in valid_masks:
        dp[mask] = 1

    # Step 3: Iterate from the second row to the n-th row.
    for i in range(1, n):
        new_dp = collections.defaultdict(int)
        # For each valid mask for the current row...
        for mask in valid_masks:
            # ...sum up the counts from compatible masks in the previous row.
            for prev_mask in valid_masks:
                # Two masks are compatible if they don't have 1s in the same column.
                if (mask & prev_mask) == 0:
                    new_dp[mask] += dp[prev_mask]
        dp = new_dp

    # Step 4: The total number of assignments is the sum of counts for all possible
    # configurations in the last row.
    total_assignments = sum(dp.values())
    
    # Step 5: Print the final result and the equation as requested.
    print(f"The problem is equivalent to counting the number of ways to place 0s and 1s on a {n}x{m} grid such that no two 1s are adjacent.")
    print("We can solve this using dynamic programming, calculating the valid assignments row by row.")
    print("\nThe number of valid assignments for the full grid is the sum of possibilities for each valid configuration of the last row:")
    
    # Get the counts for each valid mask in the final row for the equation.
    # Sorting by mask for a consistent output order.
    final_counts = [dp[mask] for mask in sorted(valid_masks)]
    
    equation_parts = [str(count) for count in final_counts]
    equation_str = " + ".join(equation_parts)
    
    print(f"\n{equation_str} = {total_assignments}")

solve_grid_assignments()