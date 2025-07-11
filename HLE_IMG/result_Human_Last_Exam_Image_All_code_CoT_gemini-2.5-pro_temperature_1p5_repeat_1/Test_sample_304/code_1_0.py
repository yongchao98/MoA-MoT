import math

def main():
    """
    Solves the graph identification puzzle by:
    1. Storing the properties of the given finite groups.
    2. Identifying the group corresponding to each of the 16 visualizations.
    3. Calculating the sum of the exponents for each of the four columns.
    """
    
    # Step 1: Define group properties (order and exponent)
    groups = {
        "PSL(2,4)":               {"order": 60, "exponent": 30},
        "Z4 x Z4":                {"order": 16, "exponent": 4},
        "D8":                     {"order": 16, "exponent": 8},
        "S4":                     {"order": 24, "exponent": 12},
        "A4":                     {"order": 12, "exponent": 6},
        "D3":                     {"order": 6,  "exponent": 6},
        "Z3 x Z3":                {"order": 9,  "exponent": 3},
        "Z2^3":                   {"order": 8,  "exponent": 2}
    }

    # Step 2: Assign groups to graphs based on analysis
    # This mapping is derived from vertex counts and structural properties of the graphs
    graph_assignments = {
        1:  "PSL(2,4)",    # Cayley, 60 vertices
        2:  "Z4 x Z4",     # Cayley, 16 vertices, abelian
        3:  "A4",          # Cayley, 12 vertices
        4:  "S4",          # Cayley, 24 vertices
        5:  "D8",          # Cayley, 16 vertices, non-abelian
        6:  "Z3 x Z3",     # Cayley, 9 vertices
        7:  "D3",          # Cayley, 6 vertices
        8:  "Z2^3",        # Cayley, 8 vertices
        9:  "PSL(2,4)",    # Cycle, 60 vertices
        10: "Z3 x Z3",     # Cycle, 9 vertices
        11: "D3",          # Cycle, 6 vertices
        12: "Z2^3",        # Cycle, 8 vertices
        13: "Z4 x Z4",     # Cycle, 16 vertices
        14: "S4",          # Cycle, 24 vertices
        15: "D8",          # Cycle, 16 vertices
        16: "A4",          # Cycle, 12 vertices
    }
    
    # The grid is arranged in a 4x4 manner.
    # Column 1: V1, V5, V9, V13
    # Column 2: V2, V6, V10, V14
    # Column 3: V3, V7, V11, V15
    # Column 4: V4, V8, V12, V16
    columns = [
        [1, 5, 9, 13],
        [2, 6, 10, 14],
        [3, 7, 11, 15],
        [4, 8, 12, 16]
    ]
    
    column_sums = []

    # Step 3: Calculate the sum of exponents for each column
    for i, col_indices in enumerate(columns):
        col_sum = 0
        exponents_in_col = []
        
        for graph_index in col_indices:
            group_name = graph_assignments[graph_index]
            exponent = groups[group_name]["exponent"]
            exponents_in_col.append(exponent)
            col_sum += exponent
        
        column_sums.append(col_sum)
        
        # Print the detailed calculation for each sum
        equation_str = " + ".join(map(str, exponents_in_col))
        print(f"S{i+1} = {equation_str} = {col_sum}")

    print(f"\nThe four column sums are: {{{', '.join(map(str, column_sums))}}}")

if __name__ == "__main__":
    main()
