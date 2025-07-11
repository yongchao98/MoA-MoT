import numpy as np

def count_assignments():
    """
    Calculates the number of valid 0/1 assignments for a 4x4 grid
    based on the given implication rules.
    """
    n = 4
    m = 4

    # Step 1: Generate all valid single-column configurations for a column of height n.
    # A configuration is valid if it has no two adjacent 1s.
    # We can check this by seeing if the bitwise AND of the number and itself shifted right is 0.
    valid_configs = []
    for i in range(1 << n):
        if (i & (i >> 1)) == 0:
            valid_configs.append(i)
    
    k = len(valid_configs)

    # Step 2: Create the transition matrix T.
    # T[i, j] = 1 if config i can be adjacent to config j, 0 otherwise.
    # Two configurations can be adjacent if their bitwise AND is 0.
    T = np.zeros((k, k), dtype=np.int64)
    for i in range(k):
        for j in range(k):
            if (valid_configs[i] & valid_configs[j]) == 0:
                T[i, j] = 1

    # Step 3: Use dynamic programming (transfer matrix method).
    # Start with dp vector for the first column (a 4x1 grid).
    # Each valid configuration is possible, so there is 1 way for each.
    dp = np.ones(k, dtype=np.int64)

    # Iterate for the remaining m-1 columns.
    for _ in range(m - 1):
        # Update the dp vector by multiplying with the transition matrix.
        # This calculates the number of ways for the next column based on the current one.
        dp = T @ dp

    # Step 4: The final dp vector contains the number of valid 4x4 grids ending
    # in each possible valid column configuration.
    final_counts = dp
    total_assignments = np.sum(final_counts)

    # Print the final summation as requested.
    equation_parts = [str(count) for count in final_counts]
    equation_str = " + ".join(equation_parts)
    print(f"The number of valid assignments for a 4x4 grid is the sum of possibilities for each valid final column configuration:")
    print(f"{equation_str} = {total_assignments}")


count_assignments()