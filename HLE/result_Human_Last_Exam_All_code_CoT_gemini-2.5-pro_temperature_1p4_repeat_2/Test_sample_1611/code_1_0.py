def solve_grid_assignments():
    """
    Calculates the number of valid 0/1 assignments for an n x m grid
    with the given implication constraints using dynamic programming.
    """
    n = 4  # Number of rows
    m = 4  # Number of columns

    # Step 1: Find all valid row masks. A mask is valid if it has no adjacent 1s.
    valid_masks = []
    for mask in range(1 << m):
        if (mask & (mask >> 1)) == 0:
            valid_masks.append(mask)

    # Step 2: Initialize DP table.
    # dp[i][mask] stores the number of ways to fill the first i+1 rows
    # where the (i+1)-th row has the configuration 'mask'.
    dp = [[0] * (1 << m) for _ in range(n)]

    # Step 3: Base case (first row).
    # For the first row, there is 1 way for each valid mask.
    for mask in valid_masks:
        dp[0][mask] = 1

    # Step 4: Fill the DP table row by row.
    for i in range(1, n):  # For rows 2 to n
        for current_mask in valid_masks:
            count = 0
            # A configuration in the current row is valid if it can be placed
            # below a valid configuration in the previous row.
            for prev_mask in valid_masks:
                # Check for vertical constraint: no two 1s in the same column.
                if (current_mask & prev_mask) == 0:
                    count += dp[i - 1][prev_mask]
            dp[i][current_mask] = count

    # Step 5: Calculate the final result.
    # The total number of assignments is the sum of ways for the last row.
    total_assignments = sum(dp[n - 1])

    # Print the final calculation as an equation
    final_sum_terms = [str(dp[n - 1][mask]) for mask in valid_masks if dp[n - 1][mask] > 0]
    equation = " + ".join(final_sum_terms)

    print(f"The number of different 0/1 assignments for a {n}x{m} grid is the sum of all valid configurations for the final row.")
    print("Final calculation:")
    print(f"{equation} = {total_assignments}")


solve_grid_assignments()
<<<1173>>>