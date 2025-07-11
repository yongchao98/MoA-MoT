import collections

def solve_grid_assignments():
    """
    Calculates the number of valid 0/1 assignments for an n x m grid
    based on the given implication rules.
    """
    n = 4
    m = 4

    # Step 1: Generate all valid row configurations for a row of length m.
    # A row is valid if it has no adjacent 1s.
    # In binary, this means checking if `(k & (k << 1)) == 0`.
    valid_rows = []
    for i in range(1 << m):
        if (i & (i << 1)) == 0:
            valid_rows.append(i)

    # Step 2: Initialize the DP table for the first row.
    # For the first row, there is exactly one way to achieve each valid configuration.
    dp = collections.defaultdict(int)
    for row in valid_rows:
        dp[row] = 1

    # Step 3: Iterate from the second row to the n-th row.
    for i in range(1, n):
        new_dp = collections.defaultdict(int)
        # For each possible valid configuration of the current row...
        for new_row in valid_rows:
            count = 0
            # ...sum the counts from compatible configurations in the previous row.
            for old_row, num_ways in dp.items():
                # Two rows are compatible if they don't have a 1 in the same column.
                # This is checked using a bitwise AND.
                if (new_row & old_row) == 0:
                    count += num_ways
            new_dp[new_row] = count
        dp = new_dp

    # Step 4: The final answer is the sum of counts for all possible last-row configurations.
    # Get the counts for each configuration of the last row.
    final_counts = list(dp.values())
    
    # Construct the equation string.
    equation_str = " + ".join(map(str, sorted(final_counts, reverse=True)))
    
    # Calculate the total sum.
    total_assignments = sum(final_counts)

    print(f"For a {n}x{m} grid, the problem is to count binary matrices with no adjacent 1s.")
    print(f"Using dynamic programming, we build the grid row by row.")
    print(f"After {n} rows, the number of valid assignments for each possible last row configuration are:")
    print(sorted(final_counts, reverse=True))
    print("\nThe total number of assignments is the sum of these values:")
    print(f"{equation_str} = {total_assignments}")

solve_grid_assignments()
<<<994>>>