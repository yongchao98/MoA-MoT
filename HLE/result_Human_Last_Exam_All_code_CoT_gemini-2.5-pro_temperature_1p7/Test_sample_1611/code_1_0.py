def solve_grid_assignments():
    """
    Calculates the number of valid 0/1 assignments for an n x m grid
    based on the given implication rules.
    """
    n = 4
    m = 4

    # The number of possible configurations (masks) for a single row.
    num_masks = 1 << m

    # Step 1: Find all valid masks for a single row.
    # A mask is valid if it does not contain adjacent 1s.
    # We check this using the bitwise operation: (mask & (mask << 1)) == 0.
    valid_masks = []
    for mask in range(num_masks):
        if (mask & (mask << 1)) == 0:
            valid_masks.append(mask)

    # Step 2: Initialize the DP table.
    # dp[i][j] will store the number of ways to fill the first i+1 rows
    # where the (i+1)-th row has the configuration of mask j.
    dp = [[0] * num_masks for _ in range(n)]

    # Step 3: Base case for the first row (i=0).
    # Any valid single-row mask is a valid assignment for the first row.
    for mask in valid_masks:
        dp[0][mask] = 1

    # Step 4: Fill the DP table for the remaining rows (i=1 to n-1).
    for i in range(1, n):
        for current_mask in valid_masks:
            count = 0
            # A transition from a previous row's mask is valid only if there
            # are no vertical conflicts (no 1s in the same column).
            # This is checked with: (current_mask & prev_mask) == 0.
            for prev_mask in valid_masks:
                if (current_mask & prev_mask) == 0:
                    count += dp[i-1][prev_mask]
            dp[i][current_mask] = count

    # Step 5: The total number of assignments is the sum of all possibilities
    # for the last row.
    total_assignments = sum(dp[n-1])

    # Print the explanation and final calculation.
    print(f"The problem is to find the number of 0/1 assignments on a {n}x{m} grid where no two adjacent (horizontal or vertical) cells are 1.")
    print("This can be solved with dynamic programming, calculating the valid assignments row by row.")
    
    # Extract the values for the final summation.
    final_row_values = [dp[n-1][mask] for mask in valid_masks]
    calculation_str = " + ".join(map(str, final_row_values))
    
    print("\nThe number of valid assignments is the sum of possibilities for the last row:")
    # The final answer is the sum of dp[3][mask] over all valid masks.
    print(f"Final calculation: {calculation_str} = {total_assignments}")


solve_grid_assignments()
print("\n<<<1234>>>")