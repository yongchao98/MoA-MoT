import numpy as np

def solve_grid_assignments():
    """
    Calculates the number of valid 0/1 assignments for a 4x4 grid
    based on the given implication rules.
    """
    n = 4  # number of rows
    m = 4  # number of columns

    # Step 1: Find all valid row configurations for a row of length m.
    # A configuration is valid if it has no adjacent 1s.
    # This can be checked with the bitwise operation (c & (c << 1)) == 0.
    valid_rows = []
    for c in range(1 << m):
        if (c & (c << 1)) == 0:
            valid_rows.append(c)
    
    num_valid_rows = len(valid_rows)
    # print(f"Found {num_valid_rows} valid row configurations: {valid_rows}")

    # Step 2: Create the transition matrix M.
    # M[i, j] = 1 if row S[i] and row S[j] are compatible (can be adjacent).
    # Compatibility means no 1s in the same column, i.e., (S[i] & S[j]) == 0.
    M = np.zeros((num_valid_rows, num_valid_rows), dtype=np.int64)
    for i in range(num_valid_rows):
        for j in range(num_valid_rows):
            if (valid_rows[i] & valid_rows[j]) == 0:
                M[i, j] = 1

    # Step 3: Use the transfer matrix method to count assignments.
    # v_k[i] = number of ways to fill first k rows, with row k having config valid_rows[i].
    
    # For the first row, each valid configuration is possible in 1 way.
    v = np.ones(num_valid_rows, dtype=np.int64)
    print(f"Number of ways for the first row (v1):\n{v}\n")

    # For subsequent rows, v_k = M * v_{k-1}.
    # We repeat this for n-1 times.
    for i in range(1, n):
        v = M.dot(v)
        print(f"Number of ways for the first {i+1} rows, categorized by the last row (v{i+1}):\n{v}\n")
        
    # Step 4: The total number of assignments is the sum of all elements in the final vector.
    total_assignments = np.sum(v)
    
    # Format the final equation as requested.
    equation_parts = [str(val) for val in v]
    equation = " + ".join(equation_parts)
    
    print("The total number of different 0/1 assignments is the sum of the elements in the final vector v4:")
    print(f"Total = {equation} = {total_assignments}")

solve_grid_assignments()