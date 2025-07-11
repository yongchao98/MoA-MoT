def solve_grid_assignments():
    """
    Calculates the number of satisfying 0/1 assignments for an n x m grid.
    """
    n = 4
    m = 4

    # Step 1 & 2: Find all valid column configurations.
    # A configuration is valid if no two vertically adjacent cells are 1.
    # We can represent a column of height n as an n-bit integer.
    valid_configs = []
    for i in range(1 << n):
        # Check for adjacent 1s: (i & (i >> 1)) must be 0
        if (i & (i >> 1)) == 0:
            valid_configs.append(i)
    
    num_valid_configs = len(valid_configs)

    # Step 3: Create the transition matrix.
    # T[i][j] = 1 if config i and config j can be adjacent.
    T = [[0] * num_valid_configs for _ in range(num_valid_configs)]
    for i in range(num_valid_configs):
        for j in range(num_valid_configs):
            # Two columns can be adjacent if they don't have 1s in the same row.
            # (config_i & config_j) must be 0
            if (valid_configs[i] & valid_configs[j]) == 0:
                T[i][j] = 1

    # Step 4: Dynamic Programming.
    # dp[i] = number of valid grids of size n x k ending with config i.
    # For the first column, there is 1 way for each valid config.
    dp = [1] * num_valid_configs

    # Iterate from the second column to the m-th column.
    for _ in range(m - 1):
        dp_new = [0] * num_valid_configs
        for j in range(num_valid_configs): # for each new state
            s = 0
            for i in range(num_valid_configs): # from each old state
                s += dp[i] * T[i][j]
            dp_new[j] = s
        dp = dp_new

    # Step 5: Calculate the total number of assignments.
    total_assignments = sum(dp)

    # Print the results in the required format.
    equation_parts = [str(x) for x in dp]
    equation = " + ".join(equation_parts)
    print(f"The number of satisfying assignments for a 4x4 grid can be found by summing the number of valid grids ending in each possible column configuration.")
    print(f"The number of valid assignments for each of the {num_valid_configs} possible configurations of the last column are:")
    print(dp)
    print("\nThe final calculation is:")
    print(f"{equation} = {total_assignments}")

solve_grid_assignments()
<<<1289>>>