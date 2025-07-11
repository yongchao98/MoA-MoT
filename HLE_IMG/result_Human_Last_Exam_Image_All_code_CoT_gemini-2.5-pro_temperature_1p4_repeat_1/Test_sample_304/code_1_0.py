import math

def lcm(a, b):
    """Computes the least common multiple of two integers."""
    return abs(a*b) // math.gcd(a, b) if a != 0 and b != 0 else 0

def group_lcm(orders):
    """Computes the lcm of a list of orders."""
    if not orders:
        return 0
    result = orders[0]
    for i in range(1, len(orders)):
        result = lcm(result, orders[i])
    return result

def solve_group_problem():
    """
    Solves the problem by identifying groups, finding exponents, and summing them.
    """
    # Step 1: Define group properties (name, element orders)
    # The exponent is the lcm of the element orders.
    group_data = {
        "PSL(2,4)": {"orders": [1, 2, 3, 5]},
        "D8":       {"orders": [1, 2, 4, 8]},
        "A4":       {"orders": [1, 2, 3]},
        "S4":       {"orders": [1, 2, 3, 4]},
        "Z4xZ4":    {"orders": [1, 2, 4]},
        "Z2^3":     {"orders": [1, 2]},
        "D3":       {"orders": [1, 2, 3]},
        "Z3xZ3":    {"orders": [1, 3]},
    }

    # Calculate exponent for each group
    group_exponents = {name: group_lcm(data["orders"]) for name, data in group_data.items()}

    # Step 2: Map each graph visualization (V1-V16) to its corresponding group
    graph_to_group = {
        1: "PSL(2,4)", 2: "D8", 3: "A4", 4: "S4",
        5: "Z4xZ4",    6: "Z2^3", 7: "D3", 8: "Z3xZ3",
        9: "PSL(2,4)", 10: "Z3xZ3", 11: "A4", 12: "D3",
        13: "Z4xZ4",   14: "S4", 15: "Z3xZ3", 16: "D8"
    }

    # Create the 4x4 grid of exponents
    grid_exponents = [[0]*4 for _ in range(4)]
    for i in range(4): # rows
        for j in range(4): # columns
            graph_num = i * 4 + j + 1
            group_name = graph_to_group[graph_num]
            exponent = group_exponents[group_name]
            grid_exponents[i][j] = exponent
            
    # Step 3: Calculate the sum of exponents for each column
    column_sums = [0, 0, 0, 0]
    for j in range(4): # columns
        for i in range(4): # rows
            column_sums[j] += grid_exponents[i][j]

    s1, s2, s3, s4 = column_sums[0], column_sums[1], column_sums[2], column_sums[3]

    # Print the results as per instructions
    print("The exponents for the groups in each column are:")
    print(f"Column 1: {grid_exponents[0][0]}, {grid_exponents[1][0]}, {grid_exponents[2][0]}, {grid_exponents[3][0]}")
    print(f"Column 2: {grid_exponents[0][1]}, {grid_exponents[1][1]}, {grid_exponents[2][1]}, {grid_exponents[3][1]}")
    print(f"Column 3: {grid_exponents[0][2]}, {grid_exponents[1][2]}, {grid_exponents[2][2]}, {grid_exponents[3][2]}")
    print(f"Column 4: {grid_exponents[0][3]}, {grid_exponents[1][3]}, {grid_exponents[2][3]}, {grid_exponents[3][3]}")
    print("\nThe sums for each column are:")
    print(f"S1 = {grid_exponents[0][0]} + {grid_exponents[1][0]} + {grid_exponents[2][0]} + {grid_exponents[3][0]} = {s1}")
    print(f"S2 = {grid_exponents[0][1]} + {grid_exponents[1][1]} + {grid_exponents[2][1]} + {grid_exponents[3][1]} = {s2}")
    print(f"S3 = {grid_exponents[0][2]} + {grid_exponents[1][2]} + {grid_exponents[2][2]} + {grid_exponents[3][2]} = {s3}")
    print(f"S4 = {grid_exponents[0][3]} + {grid_exponents[1][3]} + {grid_exponents[2][3]} + {grid_exponents[3][3]} = {s4}")
    print(f"\nThe ordered list of column sums is: {{{s1}, {s2}, {s3}, {s4}}}")
    

solve_group_problem()
<<<68, 25, 21, 29>>>