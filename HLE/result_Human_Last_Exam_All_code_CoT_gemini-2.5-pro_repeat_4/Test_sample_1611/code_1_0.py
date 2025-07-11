import numpy as np

def solve_grid_assignments():
    """
    Calculates the number of valid 0/1 assignments for an n x m grid
    based on the given implication rules.
    """
    n = 4
    m = 4

    # Step 1: Generate all valid rows of length m.
    # A row is valid if it contains no adjacent 1s.
    # This can be checked with the bitwise operation: (i & (i >> 1)) == 0.
    valid_rows = []
    for i in range(1 << m):
        if (i & (i >> 1)) == 0:
            valid_rows.append(i)
    
    num_valid_rows = len(valid_rows)

    print(f"For a {n}x{m} grid, we first find all valid configurations for a single row of length {m}.")
    print(f"There are {num_valid_rows} such valid rows.\n")

    # Step 2: Create the transfer matrix T.
    # T[i, j] = 1 if row i is compatible with row j (i.e., they can be adjacent).
    # Compatibility means their bitwise AND is 0.
    T = np.zeros((num_valid_rows, num_valid_rows), dtype=np.int64)
    for i in range(num_valid_rows):
        for j in range(num_valid_rows):
            if (valid_rows[i] & valid_rows[j]) == 0:
                T[i, j] = 1

    # Step 3: Use the transfer matrix to count assignments.
    # We start with a 1xm grid. Any valid row is a valid assignment.
    # So, we initialize a counts vector with all 1s.
    # Using object dtype to handle potentially very large numbers, though int64 is fine here.
    counts = np.ones(num_valid_rows, dtype=object)
    
    print(f"Calculating the number of assignments for a 1x{m} grid:")
    print(f"Total assignments = {np.sum(counts)}\n")

    # Step 4: Iterate for the remaining n-1 rows.
    for i in range(1, n):
        # The new counts are calculated by T * counts.
        # Since T is symmetric, T.dot(counts) works.
        counts = T.dot(counts)
        print(f"Calculating the number of assignments for a {i+1}x{m} grid:")
        print(f"Total assignments = {np.sum(counts)}\n")
        
    # Step 5: The final answer is the sum of the last counts vector.
    final_total = np.sum(counts)
    
    print("The final calculation for the 4x4 grid is the sum of the last computed counts.")
    equation_parts = [str(c) for c in counts]
    equation = " + ".join(equation_parts)
    print(f"Final sum: {equation} = {final_total}")

solve_grid_assignments()
<<<1157>>>