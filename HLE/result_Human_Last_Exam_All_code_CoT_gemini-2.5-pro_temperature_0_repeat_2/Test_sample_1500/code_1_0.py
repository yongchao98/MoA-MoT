def calculate_b2_for_sun_orbit(n, zero_simple_root_indices):
    """
    Calculates the second Betti number b2 for a coadjoint orbit of SU(n).
    The orbit is specified by which simple roots vanish on the element lambda.

    Args:
        n (int): The 'n' in SU(n). Must be >= 2.
        zero_simple_root_indices (list): A list of 1-based indices of the simple
                                         roots {alpha_1, ..., alpha_{n-1}}
                                         that are zero on lambda. An empty list
                                         means lambda is regular.
    """
    if n < 2:
        print("Error: n must be 2 or greater for SU(n).")
        return

    # The rank of the Lie algebra su(n) is n-1.
    rank_g = n - 1

    # The simple roots of su(n) (type A_{n-1}) are linearly independent.
    # The rank of the semisimple part of the stabilizer algebra, rank([g_lambda, g_lambda]),
    # is the number of simple roots that vanish on lambda.
    # We use a set to handle potential duplicate indices in the input list.
    unique_zero_roots = set(zero_simple_root_indices)
    rank_stabilizer_semisimple = len(unique_zero_roots)

    # The second Betti number b2 is given by the formula:
    # b2 = rank(g) - rank([g_lambda, g_lambda])
    b2 = rank_g - rank_stabilizer_semisimple

    # --- Output the results ---
    print(f"--- Calculation for G = SU({n}) ---")
    if not zero_simple_root_indices:
        print("Orbit type: Regular (lambda is in the interior of the Weyl alcove)")
        print("Simple roots vanishing on lambda: None")
    else:
        print("Orbit type: Singular (lambda is on the boundary of the Weyl alcove)")
        # Format the set of roots for printing
        root_str = ", ".join([f"alpha_{i}" for i in sorted(list(unique_zero_roots))])
        print(f"Simple roots vanishing on lambda: {{{root_str}}}")

    print(f"The rank of SU({n}) is n - 1 = {n} - 1 = {rank_g}")
    print(f"The rank of the semisimple part of the stabilizer is {rank_stabilizer_semisimple}")
    print(f"Resulting second Betti number b2 = {rank_g} - {rank_stabilizer_semisimple} = {b2}")
    print(f"Is b2 equal to n - 1? {'Yes' if b2 == (n - 1) else 'No'}")
    print("-" * 35 + "\n")


# --- Main execution: Demonstrate for SU(4) ---
# For SU(4), n=4 and the rank is n-1=3.
# We test if b2 is always 3.

print("Demonstrating for G = SU(4), where n-1 = 3.\n")

# Case 1: Regular lambda. No simple roots vanish.
# This corresponds to the generic coadjoint orbit G/T.
calculate_b2_for_sun_orbit(n=4, zero_simple_root_indices=[])

# Case 2: Singular lambda. The simple root alpha_2 vanishes.
# This corresponds to an orbit G/H where H is locally S(U(1)xU(2)xU(1)).
calculate_b2_for_sun_orbit(n=4, zero_simple_root_indices=[2])

# Case 3: More singular lambda. Simple roots alpha_1 and alpha_3 vanish.
# In the A_3 Dynkin diagram (o-o-o), alpha_1 and alpha_3 are not connected.
# The stabilizer is locally S(U(2)xU(2)).
calculate_b2_for_sun_orbit(n=4, zero_simple_root_indices=[1, 3])
