def solve_chromatic_number():
    """
    Calculates the maximal chromatic number of a graph G that is the sum of three cycles of length n.
    """
    # Step 1: Define the chromatic numbers for the two cases of a cycle C_n.
    # Case A: n is odd (n >= 3). The chromatic number of an odd cycle is 3.
    chi_Cn_odd = 3
    # Case B: n is even (n >= 4). The chromatic number of an even cycle is 2.
    chi_Cn_even = 2

    # Step 2: Use the property that for a graph G = G1 + G2 + G3, χ(G) = χ(G1) + χ(G2) + χ(G3).
    # Here, G = C_n + C_n + C_n, so χ(G) = 3 * χ(C_n).

    # Calculate χ(G) for odd n
    chi_G_odd = chi_Cn_odd + chi_Cn_odd + chi_Cn_odd

    # Calculate χ(G) for even n
    chi_G_even = chi_Cn_even + chi_Cn_even + chi_Cn_even

    # Step 3: Find the maximal chromatic number by comparing the two cases.
    maximal_chi_G = max(chi_G_odd, chi_G_even)

    # Step 4: Print the explanation and the final result.
    print("The graph G is the sum of three cycles of length n.")
    print("The chromatic number χ(G) is the sum of the chromatic numbers of the three cycles.")
    print("χ(G) = χ(C_n) + χ(C_n) + χ(C_n)\n")
    
    print("We have two cases for χ(C_n):")
    print(f"1. If n is odd, χ(C_n) = {chi_Cn_odd}. Then, χ(G) = {chi_Cn_odd} + {chi_Cn_odd} + {chi_Cn_odd} = {chi_G_odd}.")
    print(f"2. If n is even, χ(C_n) = {chi_Cn_even}. Then, χ(G) = {chi_Cn_even} + {chi_Cn_even} + {chi_Cn_even} = {chi_G_even}.\n")

    print(f"The maximal chromatic number is the maximum of {chi_G_odd} and {chi_G_even}, which is {maximal_chi_G}.")
    
    # As requested, output the final equation for the maximal case.
    print("\nThe final equation for the maximal case is:")
    print(f"{chi_Cn_odd} + {chi_Cn_odd} + {chi_Cn_odd} = {maximal_chi_G}")

solve_chromatic_number()