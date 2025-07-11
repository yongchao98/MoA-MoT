import collections

def solve_grid_assignments():
    """
    Calculates the number of valid 0/1 assignments for a grid
    where no two adjacent cells can be 1.
    """
    # Grid dimensions
    n = 4
    m = 4

    # Step 1: Find all valid column states for a column of height n.
    # A column state is "valid" if no two vertically adjacent cells are both 1.
    # We can represent states as integers and use bitwise operations.
    valid_states = []
    for i in range(1 << n):
        # Check if the binary representation of i has adjacent 1s.
        # (i & (i << 1)) will be non-zero if there are adjacent 1s.
        if (i & (i << 1)) == 0:
            valid_states.append(i)

    num_valid_states = len(valid_states)

    # Step 2: Initialize counts for the first column (n x 1 grid).
    # For a n x 1 grid, there is exactly one way to achieve each valid state.
    # `counts[j]` will store the number of ways a valid grid can end with `valid_states[j]`.
    counts = [1] * num_valid_states
    print(f"Problem: Counting valid 0/1 assignments on a {n}x{m} grid where no two adjacent cells are 1.")
    print(f"Number of valid states for a single column of height {n}: {num_valid_states}")
    print(f"Total assignments for a {n}x1 grid: {sum(counts)}\n")

    # Step 3: Pre-calculate compatibility between states for efficiency.
    # A transition from state_prev to state_curr is valid if (state_prev & state_curr) == 0.
    compatible_prev_indices = collections.defaultdict(list)
    for i, state_curr in enumerate(valid_states):
        for j, state_prev in enumerate(valid_states):
            if (state_curr & state_prev) == 0:
                compatible_prev_indices[i].append(j)

    # Step 4: Iterate for columns 2 to m, updating counts using dynamic programming.
    for col in range(2, m + 1):
        new_counts = [0] * num_valid_states
        for i in range(num_valid_states):
            # The number of ways for the current column `col` to be in `valid_states[i]`
            # is the sum of the counts of all compatible previous column states.
            for prev_idx in compatible_prev_indices[i]:
                new_counts[i] += counts[prev_idx]
        
        counts = new_counts
        total_for_col = sum(counts)
        print(f"Total assignments for a {n}x{col} grid: {total_for_col}")

    # Step 5: The final answer is the sum of counts for the last column.
    final_answer = sum(counts)

    print(f"\nFor the final {n}x{m} grid, the counts for each possible last column state are:")
    print(counts)
    
    final_sum_str = " + ".join(map(str, counts))
    print(f"\nThe final calculation is the sum of these counts:")
    print(f"{final_sum_str} = {final_answer}")

solve_grid_assignments()
<<<1234>>>