import math

def main():
    """
    Solves the problem by identifying groups, finding their exponents,
    and summing the exponents for each column in the grid.
    """

    # Step 1 & 2: Define the groups and their exponents.
    # The exponent is the lcm of the orders of all elements in the group.
    group_exponents = {
        "PSL(2,4)": 30,         # lcm(1,2,3,5), as PSL(2,4) is isomorphic to A5
        "Z4xZ4": 4,             # lcm(exp(Z4), exp(Z4)) = lcm(4,4)
        "D8": 8,                # Dihedral group of order 16, lcm(8,2)
        "S4": 12,               # Symmetric group on 4 letters, lcm(1,2,3,4)
        "A4": 6,                # Alternating group on 4 letters, lcm(1,2,3)
        "D3": 6,                # Dihedral group of order 6, lcm(2,3)
        "Z3xZ3": 3,             # Every non-identity element has order 3
        "Z2^3": 2,              # Every non-identity element has order 2
    }

    # Step 3: Identify the group for each visualization (V1-V16)
    # based on vertex count and graph structure, and arrange them in a 4x4 grid.
    # V1-V8 are Cayley graphs, V9-V16 are Cycle graphs.
    
    # Row 1: V1, V2, V3, V4
    # Row 2: V5, V6, V7, V8
    # Row 3: V9, V10, V11, V12
    # Row 4: V13, V14, V15, V16
    grid_assignments = [
        ["PSL(2,4)", "Z4xZ4",   "A4",      "S4"   ],  # V1(60), V2(16), V3(12), V4(24)
        ["D8",       "Z2^3",    "D3",      "Z3xZ3"],  # V5(16), V6(8),  V7(6),  V8(9)
        ["PSL(2,4)", "Z3xZ3",   "A4",      "D3"   ],  # V9(60), V10(9), V11(12),V12(6)
        ["D8",       "S4",      "Z4xZ4",   "Z2^3" ]   # V13(16),V14(24),V15(16),V16(8)
    ]
    
    # Step 4: Calculate the sum of exponents for each column.
    num_cols = 4
    column_sums = [0] * num_cols
    
    print("Calculating the sum of exponents for each column:")
    
    for j in range(num_cols):
        col_exponents = []
        for i in range(len(grid_assignments)):
            group_name = grid_assignments[i][j]
            exponent = group_exponents[group_name]
            column_sums[j] += exponent
            col_exponents.append(exponent)
        
        # Print the equation for the current column sum
        equation_str = " + ".join(map(str, col_exponents))
        print(f"S{j+1} = {equation_str} = {column_sums[j]}")

    # Step 5: Report the four column sums as an ordered list.
    print("\nThe final ordered list of the four column sums is:")
    print(column_sums)
    
    # Final answer format
    print(f"\n<<<{column_sums}>>>")

if __name__ == "__main__":
    main()
