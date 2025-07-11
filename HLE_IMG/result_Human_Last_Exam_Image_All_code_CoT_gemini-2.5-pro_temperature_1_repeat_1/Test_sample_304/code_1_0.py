import collections

def solve_graph_puzzle():
    """
    Solves the puzzle by identifying groups, finding their exponents,
    and summing the exponents for each column in the grid.
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

    # Step 2: Map each graph visualization (V1 to V16) to its group
    graph_to_group = {
        1: "PSL(2,4)",  2: "D8",       3: "A4",    4: "S4",
        5: "Z4xZ4",     6: "Z3xZ3",    7: "D3",    8: "Z2^3",
        9: "PSL(2,4)",  10: "Z3xZ3",   11: "D8",   12: "D3",
        13: "Z4xZ4",    14: "S4",      15: "Z2^3", 16: "A4",
    }

    # Step 3: Create the 4x4 grid of exponents
    exponent_grid = [[0 for _ in range(4)] for _ in range(4)]
    
    graph_number = 1
    for r in range(4):
        for c in range(4):
            group_name = graph_to_group[graph_number]
            exponent = group_exponents[group_name]
            exponent_grid[r][c] = exponent
            graph_number += 1

    # Step 4: Calculate column sums
    column_sums = [0, 0, 0, 0]
    for c in range(4):
        col_sum = 0
        for r in range(4):
            col_sum += exponent_grid[r][c]
        column_sums[c] = col_sum
        
    # Step 5: Print the results clearly
    print("The 4x4 grid of exponents is:")
    for row in exponent_grid:
        print(" ".join(map(str, row)))
    print("\nCalculating the column sums:")
    
    s1_nums = [exponent_grid[r][0] for r in range(4)]
    s1 = sum(s1_nums)
    print(f"S1 = {s1_nums[0]} + {s1_nums[1]} + {s1_nums[2]} + {s1_nums[3]} = {s1}")

    s2_nums = [exponent_grid[r][1] for r in range(4)]
    s2 = sum(s2_nums)
    print(f"S2 = {s2_nums[0]} + {s2_nums[1]} + {s2_nums[2]} + {s2_nums[3]} = {s2}")

    s3_nums = [exponent_grid[r][2] for r in range(4)]
    s3 = sum(s3_nums)
    print(f"S3 = {s3_nums[0]} + {s3_nums[1]} + {s3_nums[2]} + {s3_nums[3]} = {s3}")
    
    s4_nums = [exponent_grid[r][3] for r in range(4)]
    s4 = sum(s4_nums)
    print(f"S4 = {s4_nums[0]} + {s4_nums[1]} + {s4_nums[2]} + {s4_nums[3]} = {s4}")

    final_list = [s1, s2, s3, s4]
    print(f"\nThe final ordered list of column sums is: {final_list}")
    
    # Final answer in the required format
    print(f"\n<<<{final_list}>>>")


solve_graph_puzzle()