import collections

def count_assignments():
    """
    Calculates the number of valid 0/1 assignments for an n x m grid based on the given implications.
    The implications mean that no two adjacent (horizontally or vertically) cells can both be 1.
    This is solved using dynamic programming on a row-by-row basis.
    """
    n = 4  # Number of rows
    m = 4  # Number of columns

    # Step 1: Find all valid row masks.
    # A mask is valid if it does not contain adjacent 1s (e.g., ...11...).
    # This is checked by the condition (mask & (mask << 1)) == 0.
    num_masks = 1 << m
    valid_masks = []
    for mask in range(num_masks):
        if (mask & (mask << 1)) == 0:
            valid_masks.append(mask)

    # Step 2: Initialize DP table for the first row.
    # dp[mask] will store the number of ways to fill the grid up to the current row,
    # with the current row having the configuration 'mask'.
    # For the first row, there is exactly one way for each valid mask.
    dp = collections.defaultdict(int)
    for mask in valid_masks:
        dp[mask] = 1

    # Step 3: Fill the DP table for the remaining rows (from the 2nd to the n-th row).
    for i in range(1, n):
        new_dp = collections.defaultdict(int)
        for current_mask in valid_masks:
            count = 0
            # A valid assignment for the current row can be placed on top of
            # any valid assignment of the previous row, as long as there's no vertical conflict.
            # Vertical conflict occurs if a 1 is placed on top of another 1.
            # This is checked by the condition (current_mask & prev_mask) == 0.
            for prev_mask in valid_masks:
                if (current_mask & prev_mask) == 0:
                    count += dp[prev_mask]
            new_dp[current_mask] = count
        dp = new_dp

    # Step 4: The total number of assignments is the sum of counts for all possible last-row configurations.
    # We sort the items for a consistent output order.
    final_counts = sorted(dp.items())
    
    total = sum(dp.values())

    # Format the final equation as requested.
    equation_parts = [str(count) for mask, count in final_counts]
    equation = " + ".join(equation_parts)
    
    print(f"For a {n}x{m} grid, the total number of valid assignments is the sum of possibilities for the last row.")
    print("The number of ways to complete the grid for each valid configuration of the last row are:")
    print(f"{equation} = {total}")

count_assignments()