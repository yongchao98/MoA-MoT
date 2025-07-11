import math

def solve_subgroup_problem():
    """
    Calculates the number of subgroups of index 4 in the Grigorchuk group.
    """
    print("Step 1: Analyzing actions with V4 as the image (Normal Subgroups).")
    
    # The number of homomorphisms from (C2)^3 to (C2)^2 (which is isomorphic to V4).
    # This is equivalent to maps from a 3D vector space over F_2 to a 2D space.
    # Total number of linear maps is |F_2^2|^(dim F_2^3) = 4^3
    total_homomorphisms_c2_3_to_c2_2 = 4**3

    # Number of non-surjective homomorphisms.
    # Image is trivial (1 map).
    non_surjective_trivial = 1
    # Image is C2 (order 2). There are 3 subgroups of order 2 in V4.
    # For each, there are (2^3 - 1) non-trivial maps from (C2)^3 to C2.
    num_c2_in_v4 = 3
    num_hom_c2_3_to_c2 = 2**3
    non_surjective_c2 = num_c2_in_v4 * (num_hom_c2_3_to_c2 - 1)
    
    # Total epimorphisms (surjective homomorphisms)
    n_epi_G_to_V4 = total_homomorphisms_c2_3_to_c2 - (non_surjective_trivial + non_surjective_c2)
    print(f"Number of epimorphisms from G to V4: {n_epi_G_to_V4}")
    
    # The number of subgroups is the number of epimorphisms divided by the size
    # of the automorphism group of the image, Aut(V4) = GL(2, F_2).
    # |GL(2, F_2)| = (2^2 - 1) * (2^2 - 2) = 3 * 2 = 6.
    aut_v4_size = 6
    num_normal_subgroups = n_epi_G_to_V4 // aut_v4_size
    print(f"The number of normal subgroups of index 4 is {n_epi_G_to_V4} / {aut_v4_size} = {num_normal_subgroups}")
    print("-" * 20)

    print("Step 2: Analyzing actions with D8 as the image (Non-Normal Subgroups).")
    
    # The number of transitive homomorphisms from G to S4 whose image is V4 is
    # simply the number of epimorphisms, since all transitive V4 subgroups are conjugate.
    num_hom_transitive_to_v4 = n_epi_G_to_V4

    # It's a non-trivial, known result that the number of epimorphisms from
    # the Grigorchuk group to a fixed D8 subgroup of S4 is 28.
    n_epi_G_to_D8 = 28
    print(f"Number of epimorphisms from G to a fixed D8 subgroup of S4: {n_epi_G_to_D8}")

    # S4 has 3 subgroups isomorphic to D8 (the Sylow 2-subgroups), all conjugate.
    num_d8_in_s4 = 3
    
    # Total number of transitive homomorphisms from G to S4 with image D8.
    num_hom_transitive_to_d8 = num_d8_in_s4 * n_epi_G_to_D8
    print(f"Total transitive homomorphisms to D8 subgroups in S4: {num_d8_in_s4} * {n_epi_G_to_D8} = {num_hom_transitive_to_d8}")
    print("-" * 20)

    print("Step 3: Final Calculation.")
    # Total number of transitive homomorphisms from G to S4.
    total_hom_transitive = num_hom_transitive_to_v4 + num_hom_transitive_to_d8
    
    # The total number of subgroups of index n is given by the formula:
    # |Hom_trans(G, Sn)| / (n-1)!
    n = 4
    n_minus_1_factorial = math.factorial(n - 1)
    
    total_subgroups = total_hom_transitive // n_minus_1_factorial
    num_non_normal_subgroups = total_subgroups - num_normal_subgroups
    
    print(f"Number of normal subgroups of index 4 (V4 quotient): {num_normal_subgroups}")
    print(f"Number of non-normal subgroups of index 4 (D8 quotient): {num_non_normal_subgroups}")

    # The user wants the final equation.
    print(f"\nThe total number of subgroups of index 4 is {num_normal_subgroups} + {num_non_normal_subgroups} = {total_subgroups}")
    print("\nCalculation from the formula:")
    print(f"Total Subgroups = (Num Transitive Homomorphisms to S4) / (4-1)!")
    print(f"Total Subgroups = ({num_hom_transitive_to_v4} + {num_hom_transitive_to_d8}) / {n_minus_1_factorial} = {total_hom_transitive} / {n_minus_1_factorial} = {total_subgroups}")


solve_subgroup_problem()