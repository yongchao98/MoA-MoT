import collections

def solve_graph_puzzle():
    """
    Solves the puzzle by identifying groups, finding their exponents,
    and summing the exponents in columns.
    """

    # Step 1: Define group exponents
    group_exponents = {
        "PSL(2,4)": 30,
        "Z4 x Z4": 4,
        "D8": 8,
        "S4": 12,
        "A4": 6,
        "D3": 6,
        "Z3 x Z3": 3,
        "Z2^3": 2,
    }

    # Step 2: Map graphs to identified groups
    graph_to_group = {
        'V1': 'PSL(2,4)', 'V2': 'Z4 x Z4', 'V3': 'A4',      'V4': 'S4',
        'V5': 'D8',       'V6': 'A4',      'V7': 'D3',      'V8': 'Z2^3',
        'V9': 'Z2^3',     'V10': 'Z3 x Z3','V11': 'PSL(2,4)', 'V12': 'D3',
        'V13': 'D8',      'V14': 'S4',     'V15': 'A4',      'V16': 'Z4 x Z4',
    }

    # Create a grid of exponents based on the graph identifications
    grid = [[0]*4 for _ in range(4)]
    for i in range(4):
        for j in range(4):
            graph_id = f"V{i*4 + j + 1}"
            group_name = graph_to_group[graph_id]
            grid[i][j] = group_exponents[group_name]

    # Step 3: Calculate column sums
    column_sums = [0, 0, 0, 0]
    equations = ["", "", "", ""]
    
    for j in range(4):
        col_sum = 0
        eq_parts = []
        for i in range(4):
            col_sum += grid[i][j]
            eq_parts.append(str(grid[i][j]))
        column_sums[j] = col_sum
        equations[j] = f"S{j+1} = {' + '.join(eq_parts)} = {col_sum}"

    # Print the results
    print("The exponents for each column are summed as follows:")
    for eq in equations:
        print(eq)
    
    print("\nThe four column sums as an ordered list are:")
    print(column_sums)

solve_graph_puzzle()