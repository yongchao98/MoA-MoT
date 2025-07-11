import numpy as np

def solve_grid_assignments():
    """
    Calculates the number of valid 0/1 assignments for an n x m grid
    with the given implication constraints.
    """
    n = 4
    m = 4

    # Step 1: Find all valid row configurations for a row of length m.
    # A row is valid if it has no adjacent 1s.
    valid_rows = []
    for i in range(1 << m):
        if (i & (i << 1)) == 0:
            valid_rows.append(i)
    
    num_valid_rows = len(valid_rows)

    # Step 2: Create the transition matrix T.
    # T[i, j] = 1 if row i can be placed above row j.
    # This is true if their bitwise AND is 0.
    T = np.zeros((num_valid_rows, num_valid_rows), dtype=np.int64)
    for i in range(num_valid_rows):
        for j in range(num_valid_rows):
            if (valid_rows[i] & valid_rows[j]) == 0:
                T[i, j] = 1

    # Step 3: Calculate T^(n-1). For a 4x4 grid, this is T^3.
    if n == 1:
        total_assignments = num_valid_rows
        # We don't have a final equation for n=1, so just print the number.
        print(f"Total number of assignments for a 1x4 grid is: {total_assignments}")
        return

    Tn_minus_1 = np.linalg.matrix_power(T, n - 1)

    # Step 4: Calculate the total number of assignments.
    # This is the sum of all elements in T^(n-1).
    # We can also see the breakdown of counts for each possible last row.
    # Let u be a vector of all ones. The counts for each last row are T^(n-1) * u.
    u = np.ones(num_valid_rows, dtype=np.int64)
    assignments_per_last_row = np.dot(Tn_minus_1, u)
    
    total_assignments = np.sum(assignments_per_last_row)

    # Print the final equation as requested.
    # The components are the number of valid grids ending in each possible valid row.
    equation_parts = [str(x) for x in assignments_per_last_row]
    equation_str = " + ".join(equation_parts)
    
    print(f"The total number of assignments is the sum of counts for each possible last row configuration:")
    print(f"{equation_str} = {total_assignments}")

solve_grid_assignments()
<<<1234>>>