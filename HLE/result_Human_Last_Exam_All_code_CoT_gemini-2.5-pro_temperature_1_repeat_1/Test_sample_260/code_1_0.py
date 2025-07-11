def solve_rank_of_kernel():
    """
    Calculates the rank of the kernel of the Hurewicz homomorphism for the given space Y.
    Y is the connected sum of X1, X2, and X3.
    X1: pentagon with sides identified, pi_1(X1) = Z_5
    X2: octagon with sides identified, pi_1(X2) = Z_8
    X3: real projective plane, pi_1(X3) = Z_2
    """
    d1 = 5  # Order of pi_1(X1)
    d2 = 8  # Order of pi_1(X2)
    d3 = 2  # Order of pi_1(X3)
    n = 3   # Number of spaces in the connected sum

    print("The fundamental group of Y is G = Z_5 * Z_8 * Z_2.")
    print("The kernel K of the Hurewicz map is the commutator subgroup [G, G].")
    print("We find the rank r of K using the Euler characteristic formula chi(G) = chi(K) * chi(A),")
    print("where A is the abelianization of G.")
    print(f"The orders of the component groups are |G1|={d1}, |G2|={d2}, |G3|={d3}.")
    
    # Calculate chi(G)
    # chi(G) = (1/d1 + 1/d2 + 1/d3) - (n - 1)
    # To avoid floating point issues, we use integer arithmetic.
    # chi_G_num / chi_G_den
    chi_G_num = (d2*d3 + d1*d3 + d1*d2) - (n-1)*d1*d2*d3
    chi_G_den = d1*d2*d3
    
    # Calculate chi(A)
    # A = G_ab = Z_5 + Z_8 + Z_2, so |A| = d1*d2*d3
    chi_A_den = d1*d2*d3

    # chi(K) = 1 - r
    # chi(G) = chi(K) * chi(A)  =>  chi_G_num/chi_G_den = (1 - r) * 1/chi_A_den
    # chi_G_num = 1 - r
    # r = 1 - chi_G_num
    
    rank = 1 - chi_G_num
    
    # For printing the final equation
    term1 = d2*d3
    term2 = d1*d3
    term3 = d1*d2
    term4 = (n-1)*d1*d2*d3
    
    # This is an alternative calculation based on r = 1 - (d2d3 + d1d3 + d1d2 - (n-1)d1d2d3)
    # r = 1 - d2d3 - d1d3 - d1d2 + (n-1)d1d2d3
    # r = 1 - (d1d2 + d2d3 + d3d1) + (n-1)d1d2d3
    
    calc_term1 = d1*d2 + d2*d3 + d1*d3
    calc_term2 = (n-1) * d1*d2*d3
    final_rank = 1 - calc_term1 + calc_term2

    print("\nThe rank r can be calculated from the formula:")
    print(f"r = 1 - (|G2|*|G3| + |G1|*|G3| + |G1|*|G2|) + (n-1)*|G1|*|G2|*|G3|")
    print(f"r = 1 - ({d2}*{d3} + {d1}*{d3} + {d1}*{d2}) + ({n}-1)*{d1}*{d2}*{d3}")
    print(f"r = 1 - ({term1} + {term2} + {term3}) + {n-1}*{d1*d2*d3}")
    print(f"r = 1 - {calc_term1} + {calc_term2}")
    print(f"r = {1 - calc_term1} + {calc_term2}")
    print(f"r = {final_rank}")

solve_rank_of_kernel()