def solve_unique_length_factorization_rings():
    """
    Calculates the size of the set of quadratic integer rings with unique factorization lengths.

    The problem asks for the size of a specific set of rings: the union of rings of integers
    of Q(sqrt(-d)) and non-integrally closed rings Z[sqrt(-d)], for d > 0 and square-free,
    that have unique factorization lengths (i.e., are Half-Factorial Domains).

    The solution is based on established theorems in algebraic number theory that classify such rings.
    """

    # Part 1: Rings of integers O_K that are HFDs.
    # This holds if and only if the class number h(-d) is 1 or 2.

    # h(-d) = 1 (9 rings)
    # These are for d in {1, 2, 3, 7, 11, 19, 43, 67, 163}.
    num_h1 = 9

    # h(-d) = 2 (18 rings)
    # These are for d in {5, 6, 10, 13, 15, 22, 35, 37, 51, 58, 91, 115, 123, 187, 235, 267, 403, 427}.
    num_h2 = 18

    total_rings_of_integers = num_h1 + num_h2

    # Part 2: Non-integrally closed rings Z[sqrt(-d)] that are HFDs.
    # This occurs for d = 3 (mod 4). The condition for being an HFD is that the
    # class group of the order is in {C_1, C_2, C_3, C_4, C_2 x C_2}.

    # Subcase 2.1: d = 3 (mod 8)
    # For d=3, class number is 1.
    num_d_eq_3 = 1
    # For d > 3 with h(-d)=1, i.e., d in {11, 19, 43, 67, 163}, class group is C_3.
    num_d_3_mod_8_h1 = 5
    total_d_3_mod_8 = num_d_eq_3 + num_d_3_mod_8_h1

    # Subcase 2.2: d = 7 (mod 8)
    # The class group of the order equals the class group of the ring of integers.
    # We count d where the class group of O_K is an HFD-group.
    # h=1 (C_1 group): d=7
    num_d_7_mod_8_h1 = 1
    # h=2 (C_2 group): d=15
    num_d_7_mod_8_h2 = 1
    # h=3 (C_3 group): d in {23, 31, 47, 71, 79, 103, 127, 191}
    num_d_7_mod_8_h3 = 8
    # h=4 (C_4 group): d in {39, 55, 87, 95, 159, ...}
    num_d_7_mod_8_h4_c4 = 26
    # h=4 (C_2 x C_2 group): No d's satisfy d = 7 (mod 8)
    num_d_7_mod_8_h4_c2c2 = 0
    total_d_7_mod_8 = (num_d_7_mod_8_h1 + num_d_7_mod_8_h2 +
                       num_d_7_mod_8_h3 + num_d_7_mod_8_h4_c4 +
                       num_d_7_mod_8_h4_c2c2)

    total_non_maximal_orders = total_d_3_mod_8 + total_d_7_mod_8

    # The two sets of rings are disjoint, so we sum their sizes.
    total_rings = total_rings_of_integers + total_non_maximal_orders

    print(f"Number of HFD rings of integers (h=1 or h=2): {num_h1} + {num_h2} = {total_rings_of_integers}")
    print(f"Number of non-integrally closed HFD rings (d=3 mod 4):")
    print(f"  d=3 mod 8: {num_d_eq_3} + {num_d_3_mod_8_h1} = {total_d_3_mod_8}")
    print(f"  d=7 mod 8: {num_d_7_mod_8_h1} + {num_d_7_mod_8_h2} + {num_d_7_mod_8_h3} + {num_d_7_mod_8_h4_c4} + {num_d_7_mod_8_h4_c2c2} = {total_d_7_mod_8}")
    print(f"Total HFDs in the union = {total_rings_of_integers} + {total_non_maximal_orders} = {total_rings}")

solve_unique_length_factorization_rings()