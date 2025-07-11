import math

def lcm(a, b):
    """Computes the least common multiple of two integers."""
    return abs(a * b) // math.gcd(a, b) if a != 0 and b != 0 else 0

def calculate_sums():
    """
    Identifies groups for each graph, calculates their exponents, and sums them by column.
    """
    # Step 1: Define group exponents based on group theory.
    group_exponents = {
        "PSL(2,4)": 30, # lcm(1, 2, 3, 5)
        "Z4xZ4": 4,     # lcm(1, 2, 4)
        "D8": 8,        # lcm(1, 2, 4, 8)
        "S4": 12,       # lcm(1, 2, 3, 4)
        "A4": 6,        # lcm(1, 2, 3)
        "D3": 6,        # lcm(1, 2, 3)
        "Z3xZ3": 3,     # lcm(1, 3)
        "Z2^3": 2       # lcm(1, 2)
    }

    # Step 2: Map each graph visualization to its corresponding group.
    # This mapping is derived from analyzing the number of vertices and structure of each graph.
    graph_to_group = {
        1: "PSL(2,4)", 2: "D8", 3: "A4", 4: "S4",
        5: "Z4xZ4", 6: "Z3xZ3", 7: "D3", 8: "Z2^3",
        9: "PSL(2,4)", 10: "Z3xZ3", 11: "A4", 12: "D3",
        13: "Z4xZ4", 14: "S4", 15: "D8", 16: "Z2^3"
    }

    # Create a mapping from graph number to its group's exponent.
    graph_to_exponent = {
        graph_num: group_exponents[group_name]
        for graph_num, group_name in graph_to_group.items()
    }

    # The 4x4 grid of graph visualizations.
    grid = [
        [1, 2, 3, 4],
        [5, 6, 7, 8],
        [9, 10, 11, 12],
        [13, 14, 15, 16]
    ]

    # Step 3: Calculate the sum of exponents for each column.
    column_sums = []
    print("Calculating the sum of exponents for each column:")
    for col_index in range(4):
        s_sum = 0
        components = []
        for row_index in range(4):
            graph_num = grid[row_index][col_index]
            exponent = graph_to_exponent[graph_num]
            s_sum += exponent
            components.append(str(exponent))
        
        # Print the equation for the current column sum, as requested.
        print(f"S{col_index + 1} = {' + '.join(components)} = {s_sum}")
        column_sums.append(s_sum)

    # Final result in the specified format.
    print("\nThe final list of four column sums is:")
    print(column_sums)

calculate_sums()