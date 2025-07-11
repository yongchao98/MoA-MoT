def solve_grid_assignments():
    """
    Calculates the number of satisfying 0/1 assignments for an n x m grid
    with the given implication rules.
    This is equivalent to finding the number of independent sets in a grid graph.
    """
    n = 4
    m = 4

    # 1. Generate all valid single-column patterns (masks)
    # A pattern is valid if no two adjacent bits are 1.
    # This corresponds to the vertical implication constraint.
    num_patterns = 1 << n
    valid_patterns = []
    for p in range(num_patterns):
        if (p & (p << 1)) == 0:
            valid_patterns.append(p)

    # 2. Initialize DP table for the first column
    # The count for any valid pattern in the first column is 1.
    dp = {p: 1 for p in valid_patterns}

    # 3. Iterate through the remaining columns (from the 2nd to the m-th)
    for _ in range(m - 1):
        new_dp = {p: 0 for p in valid_patterns}
        # For each valid pattern in the current column...
        for p_curr in valid_patterns:
            # ...look at all valid patterns in the previous column.
            for p_prev, count in dp.items():
                # Check for horizontal compatibility (no 1s in the same row).
                if (p_curr & p_prev) == 0:
                    new_dp[p_curr] += count
        dp = new_dp

    # 4. The final result is the sum of counts for all patterns in the last column.
    final_counts = list(dp.values())
    total_assignments = sum(final_counts)

    # Print the final equation as requested
    equation = " + ".join(map(str, final_counts))
    print(f"The number of valid assignments for each pattern in the final column are:")
    print(list(dp.items()))
    print("\nThe final calculation is:")
    print(f"{equation} = {total_assignments}")


solve_grid_assignments()