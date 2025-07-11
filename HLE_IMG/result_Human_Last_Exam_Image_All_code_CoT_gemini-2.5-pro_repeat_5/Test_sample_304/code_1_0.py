def solve_graph_exponent_sum():
    """
    Identifies groups for each graph, finds their exponents, and sums them by column.
    """
    # Step 1: Define group exponents.
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

    # Step 2 & 3: Map grid positions to identified groups.
    # The grid is 4x4, indexed by [row][col].
    grid_groups = [
        # Col 1       Col 2       Col 3    Col 4
        ['PSL(2,4)', 'Z4xZ4',   'A4',    'S4'   ],  # Row 0: V1, V2, V3, V4
        ['D8',       'Z2^3',    'D3',    'Z3xZ3'],  # Row 1: V5, V6, V7, V8
        ['PSL(2,4)', 'Z3xZ3',   'A4',    'D3'   ],  # Row 2: V9, V10, V11, V12
        ['D8',       'S4',      'Z2^3',  'Z4xZ4']   # Row 3: V13, V14, V15, V16
    ]

    # Step 4: Calculate column sums.
    num_cols = 4
    column_sums = [0] * num_cols

    for j in range(num_cols):
        exponents_in_col = []
        for i in range(len(grid_groups)):
            group_name = grid_groups[i][j]
            exponent = group_exponents[group_name]
            column_sums[j] += exponent
            exponents_in_col.append(str(exponent))
        
        equation = f"S{j+1} = {' + '.join(exponents_in_col)} = {column_sums[j]}"
        print(equation)

    print(f"\nThe ordered list of column sums is: {column_sums}")


solve_graph_exponent_sum()
<<<[76, 21, 20, 25]>>>