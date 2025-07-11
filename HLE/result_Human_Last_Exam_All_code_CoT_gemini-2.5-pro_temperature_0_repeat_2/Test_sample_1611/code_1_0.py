import collections

def solve_grid_assignments():
    """
    Calculates the number of valid 0/1 assignments for a grid with given constraints
    using dynamic programming.
    """
    n = 4
    m = 4

    # Step 1: Generate all valid states for a single column.
    # A state is a bitmask of length n. A state is valid if no two adjacent bits are 1.
    valid_column_states = []
    for i in range(1 << n):
        if (i & (i << 1)) == 0:
            valid_column_states.append(i)

    # Step 2: Initialize DP table for the first column.
    # dp[state] = number of ways to fill the first j columns ending with 'state'.
    dp = collections.defaultdict(int)
    for state in valid_column_states:
        dp[state] = 1

    # Step 3: Iterate through columns from the second to the last.
    for j in range(1, m):
        new_dp = collections.defaultdict(int)
        # For each possible state of the current column (s2)
        for s2 in valid_column_states:
            # Sum up the counts from compatible states of the previous column (s1)
            for s1 in valid_column_states:
                # Two columns are compatible if no 1s are horizontally adjacent.
                if (s1 & s2) == 0:
                    new_dp[s2] += dp[s1]
        dp = new_dp

    # Step 4: The total number of assignments is the sum of all counts for the last column.
    final_counts = sorted(dp.items())
    
    print(f"For a {n}x{m} grid, the number of valid assignments for each possible state of the last column are:")
    
    numbers_to_sum = []
    for state, count in final_counts:
        # The state is printed as a 4-digit binary number for clarity.
        print(f"  State {state:04b}: {count}")
        numbers_to_sum.append(str(count))

    total_assignments = sum(dp.values())
    
    print("\nThe total number of assignments is the sum of these counts:")
    equation = " + ".join(numbers_to_sum)
    print(f"{equation} = {total_assignments}")


solve_grid_assignments()