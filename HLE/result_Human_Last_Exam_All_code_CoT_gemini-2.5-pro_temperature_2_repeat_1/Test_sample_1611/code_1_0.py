import collections

def solve_grid_assignments():
    """
    Calculates the number of valid 0/1 assignments for an n x m grid.
    """
    n = 4
    m = 4

    # Step 1: Find all valid rows of length m (no adjacent 1s)
    # A number i is a valid row if its binary representation has no adjacent 1s.
    # This is true if (i & (i << 1)) == 0.
    valid_rows = [i for i in range(1 << m) if (i & (i << 1)) == 0]
    num_valid_rows = len(valid_rows)
    print(f"For a grid of width {m}, there are {num_valid_rows} possible valid rows.")
    # print("Valid rows (integer representation):", valid_rows)

    # Step 2: Create a compatibility matrix for valid rows
    # T[i][j] = 1 if row i can be followed by row j, 0 otherwise.
    # Compatibility means (row_i & row_j) == 0.
    T = collections.defaultdict(list)
    for i in range(num_valid_rows):
        for j in range(num_valid_rows):
            if (valid_rows[i] & valid_rows[j]) == 0:
                T[j].append(i) # T[j] stores indices of rows compatible with row j

    # Step 3: Dynamic programming to count solutions row by row
    # counts[s_idx] is the number of valid grids of current height ending with row s_idx.
    counts = [1] * num_valid_rows

    # The number of valid 1x4 grids is the sum of initial counts.
    total_for_1xm = sum(counts)
    print(f"Number of valid assignments for a 1x4 grid: {total_for_1xm}")

    # Iterate from the 2nd row to the n-th row
    for i in range(1, n):
        new_counts = [0] * num_valid_rows
        for j in range(num_valid_rows):
            # Sum counts of all compatible previous rows
            for prev_row_idx in T[j]:
                new_counts[j] += counts[prev_row_idx]
        counts = new_counts
        
        # Print intermediate results
        if (i + 1) < n:
            print(f"Number of valid assignments for a {i+1}x{m} grid: {sum(counts)}")

    # The final answer is the sum of all counts for the n-th row
    final_answer = sum(counts)
    print(f"The calculation for the final row leads to the total sum.")
    equation_parts = " + ".join(map(str, counts))
    print(f"Final calculation: {equation_parts} = {final_answer}")
    print("\nThus, the total number of different 0/1 assignments that satisfy all implications is:")
    print(final_answer)

solve_grid_assignments()