import collections

def solve_graph_puzzle():
    """
    Solves the puzzle by identifying groups, finding their exponents, and summing them by column.
    """

    # Step 1: Define group properties (name: exponent)
    group_exponents = {
        "PSL(2,4)": 30,
        "Z4xZ4": 4,
        "D8": 8,
        "S4": 12,
        "A4": 6,
        "D3": 6,
        "Z3xZ3": 3,
        "Z2^3": 2,
    }

    # Step 2 & 3: Map visualizations to group names based on analysis
    # V1-V8 are Cayley Graphs, V9-V16 are Cycle Graphs
    visualization_to_group = {
        # Column 1
        1: "PSL(2,4)",
        5: "D8",
        9: "PSL(2,4)",
        13: "Z4xZ4",
        # Column 2
        2: "Z4xZ4",
        6: "Z3xZ3",
        10: "Z3xZ3",
        14: "A4",
        # Column 3
        3: "A4",
        7: "D3",
        11: "S4",
        15: "D8",
        # Column 4
        4: "S4",
        8: "Z2^3",
        12: "D3",
        16: "Z2^3",
    }
    
    # Organize visualizations by column
    columns = {
        1: [1, 5, 9, 13],
        2: [2, 6, 10, 14],
        3: [3, 7, 11, 15],
        4: [4, 8, 12, 16],
    }

    column_sums = []
    
    print("Calculating column sums of exponents:\n")

    # Step 4: Calculate column sums
    for i in sorted(columns.keys()):
        col_num = i
        viz_indices = columns[col_num]
        
        exponents_in_col = []
        for viz_index in viz_indices:
            group_name = visualization_to_group[viz_index]
            exponent = group_exponents[group_name]
            exponents_in_col.append(exponent)
            
        total_sum = sum(exponents_in_col)
        column_sums.append(total_sum)
        
        # Print the detailed calculation for each column
        calculation_str = " + ".join(map(str, exponents_in_col))
        print(f"S_{col_num} = {calculation_str} = {total_sum}")
        
    print("\nFinal ordered list of column sums:")
    print(column_sums)

solve_graph_puzzle()
<<<[72, 16, 32, 22]>>>