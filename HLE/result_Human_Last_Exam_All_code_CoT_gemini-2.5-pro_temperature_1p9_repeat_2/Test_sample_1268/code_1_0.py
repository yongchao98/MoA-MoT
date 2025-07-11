import math

def display_upper_bound_formula():
    """
    This script presents the formula for the upper bound of the maximum cusp height
    in relation to the covolume for Bianchi groups.

    The script assumes the user's notation (k_k,inf) refers to the height
    of the principal cusp (h_inf) for the Bianchi group associated with the
    squarefree integer N.
    """

    # --- Formula Derivation Explanation ---
    # The height of the principal cusp (h_inf) is given by:
    # h_inf = sqrt(|D_N|) / (2 * w_N)
    #
    # The covolume (V) is given by:
    # V = |D_N|^(3/2) * zeta_k(2) / (4 * pi^2)
    #
    # where zeta_k(2) = zeta(2) * L(2, chi_D_N) = (pi^2 / 6) * L(2, chi_D_N).
    #
    # By expressing |D_N| in terms of V from the second formula and substituting
    # into the first, we arrive at the relationship.
    # From the volume formula: |D_N| = ( (4 * pi^2 * V) / zeta_k(2) )^(2/3)
    # Substituting zeta_k(2): |D_N| = ( (24 * V) / L(2, chi_D_N) )^(2/3)
    #
    # Then h_inf = (1 / (2*w_N)) * sqrt(|D_N|)
    # h_inf = (1 / (2*w_N)) * ( (24 * V) / L(2, chi_D_N) )^(1/3)

    print("Based on the interpretation of (k_k,inf) as the maximal cusp height (h_inf),")
    print("the relationship with the covolume (V) for a squarefree integer N is:")
    print("\n" + "="*80)

    # --- Printing the Final Formula ---
    term1_num = 24
    term2_num = 1
    term2_den = 3
    term3_num = 2

    print(f"\n              (    {term1_num} * V    )^({term2_num}/{term2_den})")
    print( "k_k,inf  <=   ----------------------------")
    print(f"            ( {term3_num} * w_N * [L(2, chi_D_N)]^({term2_num}/{term2_den}) )")

    print("\n" + "="*80)

    # --- Explanation of Terms ---
    print("\nWhere the terms are defined as follows:")
    print("  V: The covolume of the Bianchi group PSL_2(O_N).")
    print("  N: A squarefree natural number defining the imaginary quadratic field k = Q(sqrt(-N)).")
    print("  D_N: The discriminant of the field k. It is -N if N = 3 (mod 4), and -4N otherwise.")
    print("  w_N: The number of roots of unity in the ring of integers O_N.")
    print("       - w_N = 6 if N = 3")
    print("       - w_N = 4 if N = 1")
    print("       - w_N = 2 for all other squarefree N > 1.")
    print("  L(2, chi_D_N): The value of the Dirichlet L-function for the character chi_D_N at s=2.")
    print("                 This value depends on N, but it is bounded and tends to 1 for large N.")

    # The relationship is an equality for h_inf. It serves as an upper bound as h_inf is the maximal cusp height.

if __name__ == '__main__':
    display_upper_bound_formula()
