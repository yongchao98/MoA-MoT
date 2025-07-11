import collections

def solve():
    """
    Calculates the number of satisfying 0/1 assignments for a 4x4 grid.
    """
    n = 4
    m = 4

    # 1. Find all valid row patterns for a row of length m
    valid_rows = []
    for r in range(1 << m):
        # A row is valid if it has no adjacent 1s.
        if (r & (r >> 1)) == 0:
            valid_rows.append(r)

    # dp[r] will store the number of valid grids of a certain height
    # ending with row pattern r.
    # 2. Base case: for a 1xm grid, there's 1 way for each valid row pattern.
    dp = collections.defaultdict(int)
    for r in valid_rows:
        dp[r] = 1

    # 3. DP recurrence for rows 2 to n
    for i in range(1, n):
        new_dp = collections.defaultdict(int)
        for r_curr in valid_rows:
            count = 0
            # Sum up the counts from compatible previous rows
            for r_prev, prev_count in dp.items():
                # Two rows are compatible if they have no 1s in the same column.
                if (r_curr & r_prev) == 0:
                    count += prev_count
            new_dp[r_curr] = count
        dp = new_dp

    # 4. The total number of assignments is the sum of counts for all possible last rows.
    total_assignments = sum(dp.values())

    print(f"For an {n}x{m} grid, the problem is to count binary matrices with no adjacent 1s.")
    print("\nUsing dynamic programming, we build the grid row by row.")
    print(f"After {n} rows, the number of valid assignments for each possible pattern of the last row are:")
    
    final_counts_str = []
    for r in sorted(dp.keys()):
        binary_r = format(r, f'0{m}b')
        print(f"- Row pattern {binary_r} ({r}): {dp[r]} ways")
        final_counts_str.append(str(dp[r]))
    
    print("\nThe final answer is the sum of these numbers:")
    print(" + ".join(final_counts_str) + f" = {total_assignments}")

solve()
<<<1234>>>