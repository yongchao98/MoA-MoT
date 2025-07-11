def solve():
    """
    This script calculates the number of subgroups of index 4 in the Grigorchuk group.
    The method is based on counting homomorphisms to transitive subgroups of S4.
    """

    print("The number of subgroups of index 4 is the sum of subgroups whose quotients are V4 and D4.")
    print("-" * 70)

    # Part 1: Subgroups with quotient V4 (Klein-4 group)
    print("Part 1: Calculating subgroups whose quotient is the Klein-4 group (V4)")

    # V4 is abelian, so homomorphisms G -> V4 factor through the abelianization G_ab = (Z/2Z)^3.
    # We count surjective homomorphisms (epimorphisms) from (Z/2Z)^3 to V4 = (Z/2Z)^2.
    # The number of homomorphisms from (Z/2Z)^m to (Z/2Z)^n is (2^n)^m.
    # For m=3, n=2, we have (2^2)^3 = 4^3 = 64.
    total_hom_g_ab_to_v4 = 4**3
    
    # Count non-surjective homomorphisms by considering proper subgroups of V4.
    # V4 has 3 subgroups of order 2 (isomorphic to C2).
    # 1. Homomorphisms to the trivial subgroup {e}: 1
    hom_to_trivial = 1
    # 2. Homomorphisms to one of the 3 subgroups of order 2.
    # Number of homomorphisms from (Z/2Z)^3 to a C2 subgroup: 2^3 = 8.
    # Number of surjective homomorphisms to C2: 8 - 1 (trivial map) = 7.
    # Total homomorphisms whose image is a C2 subgroup: 3 subgroups * 7 = 21.
    hom_to_c2_subgroups = 3 * (2**3 - 1)
    
    # The number of epimorphisms is the total minus the non-surjective ones.
    num_epi_g_to_v4 = total_hom_g_ab_to_v4 - hom_to_trivial - hom_to_c2_subgroups
    print(f"Number of surjective homomorphisms from G to V4: {total_hom_g_ab_to_v4} - {hom_to_trivial} - {hom_to_c2_subgroups} = {num_epi_g_to_v4}")

    # The order of the automorphism group of V4, Aut(V4), is |S3| = 6.
    aut_v4_order = 6
    print(f"Order of the automorphism group of V4: {aut_v4_order}")

    num_v4_subgroups = num_epi_g_to_v4 // aut_v4_order
    print(f"Number of subgroups of index 4 with quotient V4: {num_epi_g_to_v4} / {aut_v4_order} = {num_v4_subgroups}")
    print("-" * 70)

    # Part 2: Subgroups with quotient D4 (Dihedral group of order 8)
    print("Part 2: Calculating subgroups whose quotient is the Dihedral group (D4)")
    
    # We use the property |Hom(G, D4)| = |Hom(C2 * (C2 x C2), D4)|.
    # A homomorphism from C2 * (C2 x C2) is defined by the images of generators x, y, z.
    # The images X, Y, Z must be involutions in D4, and Y and Z must commute.
    # D4 has 6 involutions: {e, r^2, s, sr, sr^2, sr^3}.
    involutions_in_d4 = 6
    # There are 28 ordered pairs of commuting involutions (Y,Z) in D4.
    commuting_involution_pairs_d4 = 28
    
    total_hom_g_to_d4 = involutions_in_d4 * commuting_involution_pairs_d4
    print(f"Total homomorphisms from G to D4: {involutions_in_d4} * {commuting_involution_pairs_d4} = {total_hom_g_to_d4}")

    # Count non-surjective homomorphisms (image is a proper subgroup of D4).
    # Proper subgroups of D4 are 5 of type C2, 1 of type C4, and 2 of type V4.
    # The image cannot be C4.
    
    # 1. Image is trivial: 1
    hom_to_trivial_d4 = 1
    # 2. Image is C2: There are 5 C2 subgroups in D4.
    # Number of epimorphisms from G to C2 is 7. Total = 5 * 7 = 35.
    num_c2_subgroups_d4 = 5
    epi_g_to_c2 = 7
    hom_to_c2_subgroups_d4 = num_c2_subgroups_d4 * epi_g_to_c2

    # 3. Image is V4: There are 2 V4 subgroups in D4.
    # Number of epimorphisms from G to V4 is 42 (from Part 1). Total = 2 * 42 = 84.
    num_v4_subgroups_in_d4 = 2
    hom_to_v4_subgroups_d4 = num_v4_subgroups_in_d4 * num_epi_g_to_v4
    
    num_epi_g_to_d4 = total_hom_g_to_d4 - hom_to_trivial_d4 - hom_to_c2_subgroups_d4 - hom_to_v4_subgroups_d4
    print(f"Number of surjective homomorphisms from G to D4: {total_hom_g_to_d4} - ({hom_to_trivial_d4} + {hom_to_c2_subgroups_d4} + {hom_to_v4_subgroups_d4}) = {num_epi_g_to_d4}")

    # The order of the automorphism group of D4 is 8.
    aut_d4_order = 8
    print(f"Order of the automorphism group of D4: {aut_d4_order}")
    
    num_d4_subgroups = num_epi_g_to_d4 // aut_d4_order
    print(f"Number of subgroups of index 4 with quotient D4: {num_epi_g_to_d4} / {aut_d4_order} = {num_d4_subgroups}")
    print("-" * 70)

    # Part 3: Total
    print("Part 3: Total number of subgroups of index 4")
    total_subgroups = num_v4_subgroups + num_d4_subgroups
    print(f"The total number is the sum of the numbers from Part 1 and Part 2.")
    print(f"{num_v4_subgroups} + {num_d4_subgroups} = {total_subgroups}")

solve()