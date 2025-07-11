def solve_grid_assignments():
    """
    Calculates the number of valid 0/1 assignments for an n x m grid
    with the specified non-adjacency constraints for 1s.
    """
    n = 4  # Number of rows
    m = 4  # Number of columns

    num_masks = 1 << m

    # Step 1: Find all valid masks for a single row.
    # A mask is valid if it has no adjacent 1s (e.g., no '11' bit pattern).
    valid_masks = []
    for mask in range(num_masks):
        if (mask & (mask << 1)) == 0:
            valid_masks.append(mask)

    # Step 2: Initialize DP table for the first row.
    # dp[mask] stores the number of ways to tile up to the current row,
    # ending with the configuration 'mask'.
    dp = [0] * num_masks
    for mask in valid_masks:
        dp[mask] = 1

    # Step 3: Fill the DP table for the remaining rows (from row 2 to n).
    for i in range(1, n):
        new_dp = [0] * num_masks
        for current_mask in valid_masks:
            count = 0
            # Sum up the counts from compatible previous row masks.
            for prev_mask in valid_masks:
                # Compatibility check: no 1s in the same column.
                if (current_mask & prev_mask) == 0:
                    count += dp[prev_mask]
            new_dp[current_mask] = count
        dp = new_dp

    # Step 4: The final answer is the sum of all counts in the last DP row.
    total_assignments = sum(dp)
    
    # Generate the string for the final summation equation.
    final_sum_parts = []
    for mask in valid_masks:
        if dp[mask] > 0:
            final_sum_parts.append(str(dp[mask]))
    
    final_equation = " + ".join(final_sum_parts)
    
    print("The final count is the sum of valid configurations for the last row.")
    print(f"The number of assignments for each valid final row mask are: {final_sum_parts}")
    print("\nThe final equation is:")
    print(f"Total assignments = {final_equation}")
    print(f"Total assignments = {total_assignments}")

solve_grid_assignments()
<<<1234>>>