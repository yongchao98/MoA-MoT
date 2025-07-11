import collections

def solve_graph_exponent_sum():
    """
    Solves the problem by identifying groups for each graph, retrieving their exponents,
    and summing the exponents for each column in the grid.
    """
    
    # Step 1: Define group exponents
    exponents = {
        'PSL(2,4)': 30,
        'Z4xZ4': 4,
        'D8': 8,
        'S4': 12,
        'A4': 6,
        'D3': 6,
        'Z3xZ3': 3,
        'Z2^3': 2
    }

    # Step 2: Identify the group for each visualization in the 4x4 grid
    # This identification is based on vertex count and graph structure analysis.
    grid_identifications = [
        # Column 1       # Column 2     # Column 3     # Column 4
        ['PSL(2,4)',     'Z4xZ4',       'A4',          'S4'    ],  # Row 1 (V1-V4)
        ['Z2^3',         'Z3xZ3',       'D3',          'D8'    ],  # Row 2 (V5-V8)
        ['S4',           'Z3xZ3',       'A4',          'D3'    ],  # Row 3 (V9-V12)
        ['PSL(2,4)',     'Z4xZ4',       'D8',          'Z2^3'  ]   # Row 4 (V13-V16)
    ]

    # Step 3: Compute the sum of exponents for each column
    column_sums = [0, 0, 0, 0]
    
    # Store calculations for printing
    column_calculations = collections.defaultdict(list)
    
    for col_idx in range(4):
        for row_idx in range(4):
            group_name = grid_identifications[row_idx][col_idx]
            exponent = exponents[group_name]
            column_sums[col_idx] += exponent
            column_calculations[col_idx].append(exponent)
            
    # Print the detailed calculations for each column sum
    print("Calculations for the column sums:")
    for i in range(4):
        equation = ' + '.join(map(str, column_calculations[i]))
        print(f"S{i+1} = {equation} = {column_sums[i]}")
        
    # Print the final result in the specified format
    print("\nFinal ordered list of column sums:")
    print(column_sums)

solve_graph_exponent_sum()