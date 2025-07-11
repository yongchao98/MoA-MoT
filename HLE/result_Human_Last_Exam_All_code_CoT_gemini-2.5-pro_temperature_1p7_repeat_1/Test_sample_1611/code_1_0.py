import collections

def solve_grid_assignments():
    """
    Calculates the number of valid 0/1 assignments for an n x m grid
    based on the given implication rules.
    """
    n = 4
    m = 4

    # Step 1: Find all valid single-column profiles.
    # A profile is a bitmask representing the column's 0/1 values.
    # A profile is valid if it has no two vertically adjacent 1s.
    # This can be checked with the bitwise operation (p & (p >> 1)) == 0.
    valid_profiles = []
    for p in range(1 << n):
        if (p & (p >> 1)) == 0:
            valid_profiles.append(p)

    # Step 2: Initialize the DP table for the first column.
    # dp[profile] stores the number of ways to fill the grid up to the
    # current column, with the current column having that profile.
    dp = collections.defaultdict(int)
    for p in valid_profiles:
        dp[p] = 1

    # Step 3: Iterate through the remaining columns (from the second to the last).
    for j in range(1, m):
        new_dp = collections.defaultdict(int)
        for current_profile in valid_profiles:
            count = 0
            # A profile is compatible with a previous one if there are no
            # horizontally adjacent 1s. This is checked with (p1 & p2) == 0.
            for prev_profile in valid_profiles:
                if (current_profile & prev_profile) == 0:
                    count += dp[prev_profile]
            new_dp[current_profile] = count
        dp = new_dp

    # Step 4: The total number of assignments is the sum of counts for all
    # possible profiles in the last column.
    # The result needs to be sorted by profile for consistent output.
    final_counts = [dp[p] for p in sorted(dp.keys())]
    total_assignments = sum(final_counts)

    # Print the final equation as requested.
    equation_str = ' + '.join(map(str, final_counts))
    print(f"The number of ways for each possible configuration of the last column are: {sorted(dp.keys())}")
    print(f"The counts for each configuration are: {final_counts}")
    print(f"The total number of satisfying assignments is the sum of these counts:")
    print(f"{equation_str} = {total_assignments}")


solve_grid_assignments()