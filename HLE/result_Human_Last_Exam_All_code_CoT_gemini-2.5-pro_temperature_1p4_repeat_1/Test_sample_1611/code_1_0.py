import collections

def solve_grid_assignments():
    """
    Calculates the number of valid 0/1 assignments for a 4x4 grid.

    The constraints are:
    - x_i,j -> ~x_i,j+1 (no two horizontal 1s)
    - x_i,j -> ~x_i+1,j (no two vertical 1s)

    This is solved using dynamic programming row by row.
    """
    n = 4  # Number of rows
    m = 4  # Number of columns

    # 1. Generate all valid states for a single row.
    # A state is a bitmask of length m. A state is valid if no two adjacent bits are 1.
    valid_row_states = []
    for state in range(1 << m):
        if (state & (state << 1)) == 0:
            valid_row_states.append(state)

    # 2. Initialize DP for the first row.
    # dp[state] = number of ways to fill the grid up to the current row,
    # with the current row having the given state.
    dp = {state: 1 for state in valid_row_states}

    # 3. Iterate through the remaining rows (from the second row to the n-th row).
    for i in range(1, n):
        new_dp = collections.defaultdict(int)
        # For each possible state of the current row...
        for current_state in valid_row_states:
            # ...calculate the number of ways to reach it from the previous row.
            count = 0
            # A previous row state is compatible if it doesn't have a 1 in the same column.
            for prev_state, prev_count in dp.items():
                if (current_state & prev_state) == 0:
                    count += prev_count
            new_dp[current_state] = count
        dp = new_dp

    # 4. The total number of assignments is the sum of ways for all possible states in the last row.
    final_counts = list(dp.values())
    total_assignments = sum(final_counts)

    # Sort the final counts for a consistent output order
    final_counts.sort(reverse=True)
    
    # Present the result as a summation.
    equation_str = " + ".join(map(str, final_counts))
    print(f"For a {n}x{m} grid, the number of valid assignments for the final row, partitioned by the state of that row, are:")
    print(sorted(list(dp.values()), reverse=True))
    print("\nThe total number of different assignments is the sum of these values:")
    print(f"{equation_str} = {total_assignments}")


solve_grid_assignments()
<<<1234>>>