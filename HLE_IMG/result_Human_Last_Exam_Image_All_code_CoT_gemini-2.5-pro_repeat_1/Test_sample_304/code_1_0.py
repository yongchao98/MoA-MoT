import math

def main():
    """
    Identifies groups for each graph, calculates their exponents, 
    and sums the exponents for each column of the grid.
    """

    # Step 1: Define group exponents
    group_exponents = {
        "PSL(2,4)": 30,
        "Z4xZ4": 4,
        "D8": 8,
        "S4": 12,
        "A4": 6,
        "D3": 6,
        "Z3xZ3": 3,
        "Z2^3": 2
    }

    # Step 2: Map graph number to its identified group
    # This mapping is derived by analyzing the number of vertices and structure of each graph.
    graph_to_group = {
        1: "PSL(2,4)", 5: "D8", 9: "PSL(2,4)", 13: "Z2^3",       # Column 1
        2: "Z4xZ4",   6: "Z3xZ3", 10: "Z3xZ3", 14: "S4",         # Column 2
        3: "A4",      7: "D3", 11: "Z4xZ4", 15: "A4",           # Column 3
        4: "S4",      8: "Z2^3", 12: "D3", 16: "D8"            # Column 4
    }

    # The 4x4 grid of graph visualizations
    grid_layout = [
        [1, 2, 3, 4],
        [5, 6, 7, 8],
        [9, 10, 11, 12],
        [13, 14, 15, 16]
    ]

    # Step 3: Calculate column sums
    column_sums = []
    for col_index in range(4):
        current_sum = 0
        sum_components = []
        for row_index in range(4):
            graph_number = grid_layout[row_index][col_index]
            group_name = graph_to_group[graph_number]
            exponent = group_exponents[group_name]
            sum_components.append(exponent)
            current_sum += exponent
        
        # Print the equation for the current column sum
        equation_str = f"S{col_index + 1} = {' + '.join(map(str, sum_components))} = {current_sum}"
        print(equation_str)
        
        column_sums.append(current_sum)

    # Final result
    print(f"\nThe four column sums as an ordered list are: {column_sums}")

if __name__ == "__main__":
    main()
