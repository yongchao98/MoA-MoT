import math

def main():
    """
    Solves the problem by identifying groups, finding their exponents,
    and summing the exponents for each column of the grid.
    """

    # Step 1 & 2: Define groups and their exponents.
    group_exponents = {
        'PSL(2,4)': 30,
        'Z4xZ4': 4,
        'D8': 8,
        'S4': 12,
        'A4': 6,
        'D3': 6,
        'Z3xZ3': 3,
        'Z2^3': 2
    }

    # Step 1: Identify the group for each graph visualization.
    # Each group appears once as a Cayley graph (V1-V8) and once as a Cycle graph (V9-V16).
    graph_to_group = {
        'V1': 'PSL(2,4)', 'V2': 'D8',    'V3': 'A4',    'V4': 'S4',
        'V5': 'Z4xZ4',    'V6': 'Z3xZ3', 'V7': 'D3',    'V8': 'Z2^3',
        'V9': 'S4',       'V10': 'Z3xZ3', 'V11': 'A4',    'V12': 'D3',
        'V13': 'Z4xZ4',   'V14': 'PSL(2,4)', 'V15': 'Z2^3', 'V16': 'D8'
    }

    # The 4x4 grid of graph visualizations
    grid_layout = [
        ['V1', 'V2', 'V3', 'V4'],
        ['V5', 'V6', 'V7', 'V8'],
        ['V9', 'V10', 'V11', 'V12'],
        ['V13', 'V14', 'V15', 'V16']
    ]

    # Step 3 & 4: Calculate column sums.
    column_sums = []
    num_cols = len(grid_layout[0])

    print("Calculating the sum of exponents for each column:")
    for col_index in range(num_cols):
        current_sum = 0
        exponent_values = []
        for row_index in range(len(grid_layout)):
            graph_id = grid_layout[row_index][col_index]
            group_name = graph_to_group[graph_id]
            exponent = group_exponents[group_name]
            exponent_values.append(exponent)
            current_sum += exponent
        
        column_sums.append(current_sum)
        
        # Print the equation for each column sum
        equation_str = " + ".join(map(str, exponent_values))
        print(f"S{col_index + 1} = {equation_str} = {current_sum}")

    # Final result
    print("\nThe ordered list of the four column sums is:")
    print(column_sums)

if __name__ == "__main__":
    main()