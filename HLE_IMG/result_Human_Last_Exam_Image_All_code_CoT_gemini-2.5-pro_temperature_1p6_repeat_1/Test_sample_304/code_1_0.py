def solve_graph_exponents():
    """
    Solves the problem by identifying groups, finding exponents, and summing them by column.
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

    # Step 2: Map graph visualizations to groups
    # V1-V8 are Cayley graphs, V9-V16 are Cycle graphs
    graph_to_group_map = {
        1: "PSL(2,4)",
        2: "Z4 x Z4",
        3: "A4",
        4: "S4",
        5: "Z2^3",
        6: "D8",
        7: "D3",
        8: "Z3 x Z3",
        9: "PSL(2,4)",
        10: "Z3 x Z3",
        11: "A4",
        12: "D3",
        13: "Z4 x Z4",
        14: "S4",
        15: "D8",
        16: "Z2^3",
    }
    
    # Get exponents for each graph
    graph_exponents = {v: group_exponents[g] for v, g in graph_to_group_map.items()}

    # Step 3: Define grid columns
    columns = {
        1: [1, 5, 9, 13],
        2: [2, 6, 10, 14],
        3: [3, 7, 11, 15],
        4: [4, 8, 12, 16],
    }

    column_sums = []
    
    print("Calculating column sums based on group exponents:")
    print("-" * 50)

    for i in sorted(columns.keys()):
        col_graphs = columns[i]
        col_exp_values = [graph_exponents[g] for g in col_graphs]
        col_sum = sum(col_exp_values)
        column_sums.append(col_sum)
        
        # Format the equation string
        equation_str = " + ".join(map(str, col_exp_values))
        print(f"S{i} (V{',V'.join(map(str,col_graphs))}): {equation_str} = {col_sum}")

    print("-" * 50)
    print(f"The ordered list of column sums is: {column_sums}")

solve_graph_exponents()