import math

def calculate_e8_positive_roots():
    """
    Calculates the number of positive roots for the E_8 root system.
    
    The E_8 root system is defined in an 8-dimensional Euclidean space. Its roots
    are known to fall into two categories:
    1. Vectors with two non-zero coordinates, which are ±1.
    2. Vectors with eight non-zero coordinates, which are ±1/2, with a sign constraint.
    """
    rank = 8
    
    # Calculate roots of Type 1: Permutations of (±1, ±1, 0, ..., 0)
    # First, choose 2 positions for the non-zero values from the 8 available dimensions.
    num_positions = math.comb(rank, 2)
    # For each choice of two positions, the non-zero values can be (1,1), (1,-1), (-1,1), or (-1,-1).
    # This gives 4 sign combinations.
    type1_roots = num_positions * 4
    
    # Calculate roots of Type 2: Vectors of the form (±1/2, ..., ±1/2)
    # with an even number of minus signs. There are 8 coordinates.
    # The total number of sign combinations would be 2^8.
    # The constraint of having an even number of minus signs reduces this by a factor of 2.
    type2_roots = 2**(rank - 1)
    
    # The total number of roots is the sum of the two types.
    total_roots = type1_roots + type2_roots
    
    # The number of positive roots is half of the total number of roots.
    positive_roots = total_roots // 2
    
    return positive_roots, type1_roots, type2_roots, total_roots, rank

def main():
    """
    Main function to execute the task and print the solution.
    """
    N, type1, type2, total, r = calculate_e8_positive_roots()
    
    print("Step 1: Recall the relevant theorem.")
    print("The number of F_q-rational maximal tori of a reductive group G over a finite field F_q is given by q^N,")
    print("where N is the number of positive roots of the group's root system (Theorem by R. Steinberg).\n")
    
    print("Step 2: Determine the number of positive roots (N) for a group of type E_8.")
    print(f"The root system of E_8 has a rank of {r}.")
    print("The total number of roots is calculated by summing the two types of roots in its standard construction:\n")

    print(f"  - Type 1 roots are the permutations of (±1, ±1, 0, ..., 0).")
    print(f"    The number of such roots is C({r}, 2) * 4 = {math.comb(r, 2)} * 4 = {type1}.\n")
    
    print(f"  - Type 2 roots are vectors of the form (±1/2, ..., ±1/2) with an even number of minus signs.")
    print(f"    The number of such roots is 2^({r}-1) = 2^7 = {type2}.\n")

    print(f"The total number of roots in the E_8 system is the sum: {type1} + {type2} = {total}.\n")
    
    print(f"Step 3: Calculate N.")
    print(f"The number of positive roots, N, is half the total number of roots.")
    print(f"N = {total} / 2 = {N}.\n")

    print("Step 4: State the final result.")
    print("The exact number of F_q-rational maximal tori of a group of type E_8 is q^N.")
    print(f"Therefore, the number is q^{{{N}}}.")

if __name__ == "__main__":
    main()