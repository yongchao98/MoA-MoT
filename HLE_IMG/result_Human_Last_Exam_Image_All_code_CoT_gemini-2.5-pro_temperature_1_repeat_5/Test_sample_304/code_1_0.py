def solve_column_sums():
    """
    Calculates the sum of exponents for each column of the graph grid.
    """
    # Step 1: Define the exponents for each group.
    group_exponents = {
        "PSL(2,4)": 30,
        "Z4xZ4":    4,
        "D8":       8,
        "S4":       12,
        "A4":       6,
        "D3":       6,
        "Z3xZ3":    3,
        "Z2^3":     2,
    }

    # Step 2: Map the identified groups to the 4x4 grid.
    grid_groups = [
        ["PSL(2,4)", "D8",    "A4",   "S4"  ],
        ["Z4xZ4",    "Z3xZ3", "D3",   "Z2^3"],
        ["PSL(2,4)", "Z3xZ3", "A4",   "D3"  ],
        ["Z4xZ4",    "D8",    "Z2^3", "S4"  ]
    ]

    # Step 3: Create the grid of exponents from the group grid.
    exponent_grid = [
        [group_exponents[name] for name in row]
        for row in grid_groups
    ]
    
    # Step 4: Calculate and print the sum for each column.
    column_sums = []
    num_cols = 4
    num_rows = 4
    
    print("The column sums are calculated as follows:")
    
    for j in range(num_cols):
        col_exponents = [exponent_grid[i][j] for i in range(num_rows)]
        current_sum = sum(col_exponents)
        column_sums.append(current_sum)
        
        # Build and print the equation for the current column sum
        equation = " + ".join(map(str, col_exponents))
        print(f"S{j+1} = {equation} = {current_sum}")
        
    print("\nThe final ordered list of the four column sums is:")
    print(column_sums)

# Run the solver
solve_column_sums()