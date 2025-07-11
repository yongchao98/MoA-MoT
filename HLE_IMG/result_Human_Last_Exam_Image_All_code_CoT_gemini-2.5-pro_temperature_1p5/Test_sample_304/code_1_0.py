import math

def main():
    """
    Solves the puzzle by identifying groups, finding their exponents,
    and summing the exponents in columns.
    """

    # Step 1 & 2: Define group properties (exponent) and map graphs to groups.
    # The exponent is the lcm of the orders of all elements in the group.
    group_exponents = {
        "PSL(2,4)": 30,      # lcm(1, 2, 3, 5)
        "Z_4 x Z_4": 4,      # lcm(1, 2, 4)
        "D_8": 8,            # lcm(1, 2, 4, 8)
        "S_4": 12,           # lcm(1, 2, 3, 4)
        "A_4": 6,            # lcm(1, 2, 3)
        "D_3": 6,            # lcm(1, 2, 3)
        "Z_3 x Z_3": 3,      # lcm(1, 3)
        "Z_2^3": 2,          # lcm(1, 2)
    }

    # Based on visual analysis of the graphs
    graph_to_group = {
        1: "PSL(2,4)",   2: "Z_4 x Z_4",   3: "A_4",         4: "S_4",
        5: "Z_2^3",      6: "D_8",         7: "D_3",         8: "Z_3 x Z_3",
        9: "PSL(2,4)",  10: "Z_2^3",      11: "Z_4 x Z_4",  12: "D_3",
       13: "Z_3 x Z_3", 14: "S_4",        15: "A_4",        16: "D_8"
    }

    # Step 3: Construct the 4x4 grid of graph numbers
    grid_layout = [
        [1,  2,  3,  4],
        [5,  6,  7,  8],
        [9, 10, 11, 12],
        [13, 14, 15, 16]
    ]

    # Create the grid of exponents
    exponent_grid = [[0] * 4 for _ in range(4)]
    for i in range(4):
        for j in range(4):
            graph_num = grid_layout[i][j]
            group_name = graph_to_group[graph_num]
            exponent = group_exponents[group_name]
            exponent_grid[i][j] = exponent

    # Step 4: Calculate column sums
    column_sums = [0] * 4
    for j in range(4):
        sum_eq_str = f"S_{j+1} = "
        for i in range(4):
            column_sums[j] += exponent_grid[i][j]
            sum_eq_str += str(exponent_grid[i][j])
            if i < 3:
                sum_eq_str += " + "
        sum_eq_str += f" = {column_sums[j]}"
        print(sum_eq_str)
        
    print("\nThe four column sums are:")
    print(f"{{{', '.join(map(str, column_sums))}}}")
    
    # Final answer in the required format
    final_answer = f"{{{', '.join(map(str, column_sums))}}}"
    print(f"\n<<<{final_answer}>>>")

if __name__ == "__main__":
    main()