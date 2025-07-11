def solve_grid_assignments():
    """
    Calculates the number of satisfying 0/1 assignments for a 4x4 grid
    with the given implication constraints.
    """
    n = 4  # Number of rows
    m = 4  # Number of columns

    # Step 1: Find all valid row patterns for a single row of length m.
    # A row pattern is valid if it has no two adjacent 1s.
    # We can check this by seeing if the bitwise AND of the pattern and the
    # pattern shifted right by one is zero.
    valid_patterns = []
    for i in range(1 << m):
        if (i & (i >> 1)) == 0:
            valid_patterns.append(i)

    num_patterns = len(valid_patterns)

    # Step 2: Build the transition matrix T.
    # T[i][j] = 1 if pattern i and pattern j are compatible, otherwise 0.
    # Two patterns are compatible if their bitwise AND is 0 (no 1s in the same column).
    T = [[0] * num_patterns for _ in range(num_patterns)]
    for i in range(num_patterns):
        for j in range(num_patterns):
            p1 = valid_patterns[i]
            p2 = valid_patterns[j]
            if (p1 & p2) == 0:
                T[i][j] = 1

    # Step 3: Use the transfer matrix to count assignments.
    # 'counts[i]' stores the number of ways to form a valid grid of a certain
    # height ending with the i-th valid pattern.
    # For the first row, there is 1 way for each valid pattern.
    counts = [1] * num_patterns

    # Iterate for the remaining n-1 rows.
    for _ in range(n - 1):
        new_counts = [0] * num_patterns
        # new_counts = T * counts
        for i in range(num_patterns):
            for j in range(num_patterns):
                new_counts[i] += T[i][j] * counts[j]
        counts = new_counts

    # Step 4: The total number is the sum of counts for all patterns in the last row.
    total_assignments = sum(counts)
    
    # The problem asks to output the final equation.
    # The `counts` vector contains the number of valid 4x4 grids ending
    # in each of the possible valid final rows.
    equation_parts = [str(c) for c in counts]
    equation_str = " + ".join(equation_parts)
    print(f"{equation_str} = {total_assignments}")


solve_grid_assignments()