import math

def get_exponents():
    """
    Returns a dictionary mapping group names to their exponents.
    """
    return {
        "PSL(2,4)": 30,   # Projective Special Linear Group PSL(2,4)
        "Z4_Z4": 4,       # Direct Product of two Cyclic Groups of order 4
        "D_8": 8,         # Dihedral Group of order 16
        "S_4": 12,        # Symmetric Group of degree 4
        "A_4": 6,         # Alternating Group of degree 4
        "D_3": 6,         # Crystallographic Point Group D_3
        "Z3_Z3": 3,       # Direct Product of two Cyclic Groups of order 3
        "Z2_3": 2         # Elementary Abelian Group of order 8
    }

def get_grid_map():
    """
    Returns the mapping of graph visualizations (V_i) to group names.
    This mapping is derived by analyzing the number of vertices and structure of each graph.
    """
    return {
        # Col 1          # Col 2          # Col 3          # Col 4
        1: "S_4",        2: "D_8",        3: "D_3",        4: "PSL(2,4)",
        5: "Z2_3",       6: "Z3_Z3",      7: "A_4",        8: "Z4_Z4",
        9: "PSL(2,4)",   10: "Z2_3",      11: "S_4",       12: "D_3",
        13: "Z3_Z3",     14: "Z4_Z4",     15: "A_4",       16: "D_8"
    }

def main():
    """
    Calculates and prints the sum of exponents for each column of the grid.
    """
    exponents = get_exponents()
    grid_map = get_grid_map()

    # Calculate column sums
    columns = [[1, 5, 9, 13], [2, 6, 10, 14], [3, 7, 11, 15], [4, 8, 12, 16]]
    column_sums = []
    
    for i, col_indices in enumerate(columns):
        col_exponents = [exponents[grid_map[idx]] for idx in col_indices]
        col_sum = sum(col_exponents)
        column_sums.append(col_sum)
        
        # Format the output string for the equation
        equation_str = " + ".join(map(str, col_exponents))
        print(f"S{i+1} = {equation_str} = {col_sum}")

    # Print the final result in the specified format
    result_str = ", ".join(map(str, column_sums))
    print(f"The ordered list of column sums is: {{{result_str}}}")

if __name__ == "__main__":
    main()
