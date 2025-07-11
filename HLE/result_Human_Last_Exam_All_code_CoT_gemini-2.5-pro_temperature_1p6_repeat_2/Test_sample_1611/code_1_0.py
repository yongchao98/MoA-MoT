import numpy as np

def solve():
    """
    Calculates the number of valid 0/1 assignments for an n x m grid
    based on the given implications.
    """
    n = 4
    m = 4

    # Step 1: Find all valid row configurations for an m-bit row.
    # A row is valid if it does not contain adjacent 1s.
    # We can check this by seeing if the bitwise AND of the number and itself shifted left is zero.
    valid_rows = []
    for i in range(1 << m):
        if (i & (i << 1)) == 0:
            valid_rows.append(i)
    
    num_valid_rows = len(valid_rows)
    print(f"For a row of length {m}, there are {num_valid_rows} valid configurations.")
    # print(f"Valid row patterns (in decimal): {valid_rows}")
    # Example patterns: 0 (0000), 1 (0001), 2 (0010), 5 (0101), 10 (1010)

    # Step 2: Create the compatibility matrix T.
    # T[i, j] is 1 if row i and row j are compatible (can be adjacent).
    # Compatibility means they have no 1s in the same column, i.e., (row_i & row_j) == 0.
    T = np.zeros((num_valid_rows, num_valid_rows), dtype=int)
    for i in range(num_valid_rows):
        for j in range(num_valid_rows):
            if (valid_rows[i] & valid_rows[j]) == 0:
                T[i, j] = 1

    # Step 3: Use dynamic programming (matrix exponentiation) to find the result.
    # Start with a dp vector for a 1-row grid. Each valid row is one possible configuration.
    dp = np.ones(num_valid_rows, dtype=int)
    print(f"\nNumber of ways for a 1x{m} grid: {np.sum(dp)}")

    # We iterate from the 2nd row to the n-th row.
    # For each row, the new dp vector is T multiplied by the previous dp vector.
    for i in range(1, n):
        dp = T @ dp
        print(f"Number of ways for a {i+1}x{m} grid: {np.sum(dp)}")

    # Step 4: The final answer is the sum of the last dp vector.
    final_counts = dp
    total_assignments = np.sum(final_counts)
    
    print("\nThe number of ways to form a valid 4x4 grid ending with each possible valid row are:")
    print(final_counts)
    
    # Print the final sum as an equation
    final_equation = " + ".join(map(str, final_counts))
    print(f"\nThe total number of different assignments is the sum of these counts:")
    print(f"{final_equation} = {total_assignments}")

solve()
<<<1234>>>