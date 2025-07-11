import collections

def solve_grid_assignments():
    """
    Calculates the number of valid 0/1 assignments for an n x m grid
    based on the given implication rules.
    """
    # Grid dimensions
    n = 4
    m = 4

    # Step 1: Find all valid row assignments for a row of length m.
    # A row is valid if it has no two adjacent 1s.
    valid_rows = []
    for i in range(1 << m):
        # Check for adjacent 1s by shifting and ANDing.
        # If (i & (i << 1)) is 0, it means no '11' pattern exists.
        if (i & (i << 1)) == 0:
            valid_rows.append(i)
    
    # Sort for consistent ordered output
    valid_rows.sort()

    # Step 2: Dynamic programming to count assignments.
    # dp[mask] stores the number of ways to fill the first `i` rows,
    # ending with the row `mask`. We initialize for the first row (i=1).
    dp = collections.OrderedDict()
    for r in valid_rows:
        dp[r] = 1

    # Step 3: Iterate from the second row to the n-th row.
    for i in range(2, n + 1):
        new_dp = collections.OrderedDict()
        for r_curr in valid_rows:
            count = 0
            # A previous row `r_prev` is compatible with the current row `r_curr`
            # if there are no 1s in the same column, i.e., (r_curr & r_prev) == 0.
            for r_prev, prev_count in dp.items():
                if (r_curr & r_prev) == 0:
                    count += prev_count
            new_dp[r_curr] = count
        dp = new_dp

    # Step 4: The total number of assignments is the sum of counts for the last row.
    total_assignments = sum(dp.values())

    # Format the final output to show the summation.
    sum_parts = [str(v) for v in dp.values()]
    equation = " + ".join(sum_parts)
    
    print(f"For a {n}x{m} grid, the number of ways for the final row to have a specific valid configuration are:")
    for mask, count in dp.items():
        print(f"  - Row configuration {mask:0{m}b}: {count} ways")
    
    print("\nThe total number of different 0/1 assignments is the sum of these possibilities:")
    print(f"Total = {equation}")
    print(f"Total = {total_assignments}")

solve_grid_assignments()