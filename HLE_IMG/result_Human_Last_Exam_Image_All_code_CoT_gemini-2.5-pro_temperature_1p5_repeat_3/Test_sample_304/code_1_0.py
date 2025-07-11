import collections

def solve_graph_exponents():
    """
    This function identifies the group for each graph, finds its exponent,
    and calculates the sum of exponents for each column in the grid.
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

    # Step 2: Map graph visualization index to its corresponding group name
    # This mapping is based on analysis of the graph's order (number of vertices) and structure
    graph_to_group = {
        # Column 1
        1: "PSL(2,4)", 5: "Z4xZ4", 9: "PSL(2,4)", 13: "D8",
        # Column 2
        2: "D8", 6: "Z3xZ3", 10: "Z3xZ3", 14: "S4",
        # Column 3
        3: "A4", 7: "D3", 11: "A4", 15: "Z4xZ4",
        # Column 4
        4: "S4", 8: "Z2^3", 12: "D3", 16: "Z2^3"
    }

    # Retrieve exponent for each graph
    exponents = {v_idx: group_exponents[group_name] for v_idx, group_name in graph_to_group.items()}

    # Step 3: Define grid columns
    columns = {
        "S1": [1, 5, 9, 13],
        "S2": [2, 6, 10, 14],
        "S3": [3, 7, 11, 15],
        "S4": [4, 8, 12, 16],
    }

    # Step 4: Calculate and print sums for each column
    column_sums = collections.OrderedDict()
    print("Calculation of column sums:")
    for col_name, v_indices in columns.items():
        col_exponents = [exponents[i] for i in v_indices]
        col_sum = sum(col_exponents)
        column_sums[col_name] = col_sum
        
        # Format the equation string
        equation_parts = [str(exp) for exp in col_exponents]
        equation_str = " + ".join(equation_parts)
        
        print(f"{col_name} = {equation_str} = {col_sum}")
        
    final_sums_list = list(column_sums.values())
    print("\nFinal ordered list of column sums:")
    print(final_sums_list)

solve_graph_exponents()