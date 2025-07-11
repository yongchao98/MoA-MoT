import collections

def solve_grid_assignments():
    """
    Calculates the number of valid 0/1 assignments for an n x m grid.
    The rules are: no two adjacent cells (horizontally or vertically) can both be 1.
    """
    n = 4
    m = 4

    # Step 1: Find all valid row configurations (masks)
    # A row is valid if no two adjacent bits are 1.
    # This can be checked with the condition (s & (s << 1)) == 0.
    valid_rows = []
    for s in range(1 << m):
        if (s & (s << 1)) == 0:
            valid_rows.append(s)

    # Step 2: Initialize DP table for the first row
    # dp[s] stores the number of ways to fill the grid up to the current row,
    # ending with the row configuration 's'.
    dp = collections.defaultdict(int)
    for s in valid_rows:
        dp[s] = 1

    # Step 3: Iterate through the rows from the second to the last
    for i in range(1, n):
        new_dp = collections.defaultdict(int)
        # For each valid configuration of the current row
        for s_curr in valid_rows:
            # Sum up the counts from compatible previous row configurations
            for s_prev in valid_rows:
                # Check for vertical conflicts (a 1 above a 1)
                if (s_curr & s_prev) == 0:
                    new_dp[s_curr] += dp[s_prev]
        dp = new_dp

    # Step 4: The final DP table `dp` contains the counts for the n-th row.
    # The total is the sum of these counts.
    final_counts = [dp[s] for s in valid_rows]
    total_assignments = sum(final_counts)
    
    # Step 5: Print the result as a sum equation.
    equation_str = " + ".join(map(str, final_counts))
    print(f"The total number of assignments is the sum of possibilities for each valid configuration of the last row:")
    print(f"{equation_str} = {total_assignments}")

solve_grid_assignments()