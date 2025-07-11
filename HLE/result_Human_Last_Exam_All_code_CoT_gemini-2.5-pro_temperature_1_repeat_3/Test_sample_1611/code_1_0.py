import collections

def solve_grid_assignments():
    """
    Calculates the number of valid 0/1 assignments for a grid
    where no two adjacent cells can both be 1.
    """
    n = 4  # Number of rows
    m = 4  # Number of columns

    # Step 1: Find all valid row masks for a row of width m.
    # A mask is valid if no two adjacent bits are 1.
    # The condition (mask & (mask << 1)) == 0 checks this efficiently.
    valid_masks = []
    for mask in range(1 << m):
        if (mask & (mask << 1)) == 0:
            valid_masks.append(mask)

    # Step 2: Initialize DP table for the first row.
    # dp[mask] stores the number of valid assignments for the first i rows,
    # with the i-th row having the configuration 'mask'.
    dp = collections.defaultdict(int)
    for mask in valid_masks:
        dp[mask] = 1

    # Step 3: Iterate through the remaining rows (from 2 to n).
    for i in range(1, n):
        new_dp = collections.defaultdict(int)
        # For each possible valid mask for the current row (i+1)
        for new_mask in valid_masks:
            count = 0
            # Sum up the counts from the previous row (i)
            # for all compatible masks.
            # Compatibility means (new_mask & old_mask) == 0.
            for old_mask in valid_masks:
                if (new_mask & old_mask) == 0:
                    count += dp[old_mask]
            new_dp[new_mask] = count
        dp = new_dp

    # Step 4: Calculate the total number of assignments by summing up
    # the counts for all possible masks in the final row.
    total_assignments = sum(dp.values())

    print(f"For a {n}x{m} grid, the problem is to find the number of assignments where no two adjacent (horizontally or vertically) cells are 1.")
    print("This can be solved with dynamic programming, row by row.")
    print(f"The number of ways to fill a {n}x{m} grid is the sum of ways to form a valid 4th row, given the valid configurations of the first 3 rows.")

    # Sort the final DP table items by mask for a consistent output order
    final_counts = sorted(dp.items())
    sum_string = " + ".join([str(v) for k, v in final_counts])

    print("\nThe final calculation is the sum of counts for each valid configuration of the last row:")
    print(f"{sum_string} = {total_assignments}")

solve_grid_assignments()