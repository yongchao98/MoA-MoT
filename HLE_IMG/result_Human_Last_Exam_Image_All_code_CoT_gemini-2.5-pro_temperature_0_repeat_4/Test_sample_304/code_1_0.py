import collections

def solve_graph_exponent_sum():
    """
    Identifies groups for each graph, finds their exponents, and sums them by column.
    """
    # Step 1: Define the properties (order and exponent) for each group.
    groups = {
        "PSL(2,4)": {"order": 60, "exponent": 30},
        "Z4xZ4": {"order": 16, "exponent": 4},
        "D8": {"order": 16, "exponent": 8},
        "S4": {"order": 24, "exponent": 12},
        "A4": {"order": 12, "exponent": 6},
        "D3": {"order": 6, "exponent": 6},
        "Z3xZ3": {"order": 9, "exponent": 3},
        "Z2^3": {"order": 8, "exponent": 2},
    }

    # Step 2: Identify the group for each graph visualization based on its order and structure.
    # The graph number is the key, and the group name is the value.
    graph_to_group = {
        1: "PSL(2,4)", 5: "D8", 9: "S4", 13: "Z4xZ4",         # Column 1
        2: "Z4xZ4",   6: "Z2^3", 10: "Z3xZ3", 14: "PSL(2,4)", # Column 2
        3: "A4",      7: "D3", 11: "A4", 15: "D8",            # Column 3
        4: "S4",      8: "Z3xZ3", 12: "D3", 16: "Z2^3"         # Column 4
    }

    # Create a dictionary to hold the exponent for each graph
    graph_exponents = {i: groups[graph_to_group[i]]["exponent"] for i in range(1, 17) if i in graph_to_group}

    # Step 3: Define the graphs in each column.
    columns = {
        1: [1, 5, 9, 13],
        2: [2, 6, 10, 14],
        3: [3, 7, 11, 15],
        4: [4, 8, 12, 16]
    }

    # Step 4: Calculate the sum of exponents for each column.
    column_sums = collections.OrderedDict()
    for col_num, graph_ids in sorted(columns.items()):
        exponents_in_col = [graph_exponents[gid] for gid in graph_ids]
        column_sums[col_num] = sum(exponents_in_col)

    # Step 5: Report the results.
    print("Group assignments and their exponents:")
    for i in range(1, 5):
        print(f"\nColumn {i}:")
        for graph_id in columns[i]:
            group_name = graph_to_group[graph_id]
            exponent = graph_exponents[graph_id]
            print(f"  V_{graph_id}: {group_name}, Exponent = {exponent}")

    print("\nCalculating column sums:")
    s1_exps = [graph_exponents[gid] for gid in columns[1]]
    s2_exps = [graph_exponents[gid] for gid in columns[2]]
    s3_exps = [graph_exponents[gid] for gid in columns[3]]
    s4_exps = [graph_exponents[gid] for gid in columns[4]]
    
    s1 = sum(s1_exps)
    s2 = sum(s2_exps)
    s3 = sum(s3_exps)
    s4 = sum(s4_exps)

    print(f"S1 = {' + '.join(map(str, s1_exps))} = {s1}")
    print(f"S2 = {' + '.join(map(str, s2_exps))} = {s2}")
    print(f"S3 = {' + '.join(map(str, s3_exps))} = {s3}")
    print(f"S4 = {' + '.join(map(str, s4_exps))} = {s4}")

    final_sums = list(column_sums.values())
    print(f"\nThe four column sums as an ordered list are: {final_sums}")

solve_graph_exponent_sum()