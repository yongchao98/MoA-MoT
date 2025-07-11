import numpy as np

def solve_grid_assignments():
    """
    Calculates the number of valid 0/1 assignments for an n x m grid.
    """
    # Grid dimensions
    n, m = 4, 4

    # Step 1: Find all valid states for a single column.
    # A state is an n-bit integer. It's valid if it has no adjacent 1s.
    valid_states = []
    for i in range(1 << n):
        if (i & (i << 1)) == 0:
            valid_states.append(i)
    
    num_states = len(valid_states)

    # Step 2: Build the transition matrix T.
    # T[i, j] = 1 if state i can be followed by state j.
    # This is true if (state_i & state_j) == 0.
    T = np.zeros((num_states, num_states), dtype=int)
    for i in range(num_states):
        for j in range(num_states):
            if (valid_states[i] & valid_states[j]) == 0:
                T[i, j] = 1

    # Step 3: Calculate the number of assignments iteratively.
    # 'counts' is a vector where counts[i] is the number of valid grids of
    # a certain width ending with the i-th valid state.
    # For a 1-column grid, there is one way for each valid state.
    # We use object dtype to handle potentially large numbers, though int is fine here.
    counts = np.ones(num_states, dtype=object)

    # To get from a grid of width k to k+1, we multiply by T.
    # We do this m-1 times to get to a grid of width m.
    # Since we need T^T @ counts, and T is symmetric, we can just use T @ counts.
    for i in range(1, m):
        counts = T @ counts

    # The final 'counts' vector contains the number of valid 4x4 assignments
    # for each possible final column state.
    
    # Step 4: Print the final calculation.
    # The total number is the sum of all values in the final counts vector.
    total_assignments = np.sum(counts)
    
    print(f"For a {n}x{m} grid, the number of valid assignments is given by the sum:")
    equation = " + ".join(map(str, counts))
    print(f"{equation} = {total_assignments}")

solve_grid_assignments()