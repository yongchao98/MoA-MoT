import numpy as np

def solve_graph_exponent_puzzle():
    """
    Solves the puzzle by identifying groups for each graph, finding their exponents,
    and summing the exponents for each column in the grid.
    """
    
    # Step 1: Define the groups and their exponents.
    exponents = {
        'PSL(2,4)': 30,
        'Z4xZ4': 4,
        'D8': 8,
        'S4': 12,
        'A4': 6,
        'D3': 6,
        'Z3xZ3': 3,
        'Z2^3': 2,
    }

    # Step 2: Assign each graph in the 4x4 grid to a group based on analysis.
    # The grid is represented as a list of lists (rows).
    # V1-V4 form the first row, V5-V8 the second, and so on.
    grid_assignments = [
        # Col 1,     Col 2,     Col 3,   Col 4
        ['S4',    'Z4xZ4',   'D3',    'PSL(2,4)'], # Row 1: V1, V2, V3, V4
        ['Z3xZ3', 'D8',      'A4',    'Z2^3'],     # Row 2: V5, V6, V7, V8
        ['S4',    'Z3xZ3',   'A4',    'D3'],       # Row 3: V9, V10, V11, V12
        ['Z2^3',  'PSL(2,4)','Z4xZ4', 'D8']        # Row 4: V13, V14, V15, V16
    ]

    # Convert to a NumPy array for easier column slicing
    grid_data = np.array(grid_assignments)
    
    # Step 3: Calculate the sum of exponents for each column.
    column_sums = []
    print("Calculating the sum of exponents for each column:")
    
    for col_idx in range(4):
        column_groups = grid_data[:, col_idx]
        current_sum = 0
        sum_equation = []
        
        for group_name in column_groups:
            exponent_value = exponents[group_name]
            current_sum += exponent_value
            sum_equation.append(str(exponent_value))
            
        column_sums.append(current_sum)
        
        print(f"S{col_idx+1} = {' + '.join(sum_equation)} = {current_sum}")
        
    # Step 4: Report the final ordered list of sums.
    print("\nThe four column sums as an ordered list are:")
    print(column_sums)

# Run the solver
solve_graph_exponent_puzzle()