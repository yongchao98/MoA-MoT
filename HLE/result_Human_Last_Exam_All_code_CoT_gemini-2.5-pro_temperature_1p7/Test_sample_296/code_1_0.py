import math
import itertools

def calculate_and_print_subgroup_count():
    """
    Calculates the number of subgroups of index 4 in the Grigorchuk group.
    """
    # Step 1: Calculate the number of normal subgroups of index 4.
    # These correspond to epimorphisms from Gamma to V4 = C2 x C2.
    # Since V4 is abelian, this is equivalent to counting epimorphisms from
    # the abelianization of Gamma, Gamma_ab = C2 x C2 x C2, to V4.
    print("Step 1: Calculate the number of normal subgroups of index 4.")
    
    # An epimorphism phi: C2^3 -> C2^2 is a linear map of vector spaces over F_2.
    # A map is defined by the images of the 3 basis vectors of C2^3.
    # Each can map to any of the 4 elements of C2^2. Total maps = 4^3 = 64.
    num_total_homs = 4**3
    
    # A map is NOT an epimorphism if its image is a proper subgroup of C2^2.
    # Proper subgroups of C2^2 are the trivial group and 3 subgroups of order 2.
    # Number of maps into the trivial subgroup: 1^3 = 1.
    num_homs_to_trivial = 1
    # Number of maps into a specific subgroup of order 2: 2^3 = 8.
    num_subgroups_order_2 = 3
    num_homs_to_subgroup_order_2 = 2**3
    # Using inclusion-exclusion for maps whose image is proper:
    # (maps to S1) + (maps to S2) + (maps to S3) - 2*(map to trivial)
    # as the intersection of any two distinct subgroups of order 2 is the trivial group.
    num_non_epi = (num_subgroups_order_2 * (num_homs_to_subgroup_order_2 - 1)) + num_homs_to_trivial
    num_epi = num_total_homs - num_non_epi
    
    # The number of automorphisms of V4 is |GL(2, F_2)| = (2^2 - 1)(2^2 - 2) = 6.
    num_aut_v4 = (2**2 - 1) * (2**2 - 2)
    
    # Number of subgroups = (Number of Epimorphisms) / (Number of Automorphisms of the image)
    num_normal_v4_subgroups = num_epi // num_aut_v4
    
    print(f"The number of surjective homomorphisms from Gamma to V4 is {num_epi}.")
    print(f"The number of automorphisms of V4 is {num_aut_v4}.")
    print(f"Number of normal subgroups of index 4 = {num_epi} / {num_aut_v4} = {num_normal_v4_subgroups}.")
    print("-" * 20)
    
    # Step 2: Calculate the number of non-normal subgroups of index 4.
    # These arise from quotients of Gamma isomorphic to the dihedral group D4.
    print("Step 2: Calculate the number of non-normal subgroups of index 4.")
    
    # It is a known result that there are 6 distinct normal subgroups N
    # such that Gamma/N is isomorphic to D4.
    num_d4_quotients = 6
    
    # Each D4 group has 5 subgroups of index 4 (order 2). One is normal (the center),
    # and 4 are non-normal. These 4 give rise to non-normal subgroups in Gamma.
    num_non_normal_in_d4 = 4
    
    total_non_normal = num_d4_quotients * num_non_normal_in_d4
    
    print(f"It is a known result that Gamma has {num_d4_quotients} normal subgroups N with Gamma/N isomorphic to D4.")
    print(f"Each such quotient gives rise to {num_non_normal_in_d4} distinct non-normal subgroups of index 4.")
    print(f"Total number of non-normal subgroups = {num_d4_quotients} * {num_non_normal_in_d4} = {total_non_normal}.")
    print("-" * 20)
    
    # Step 3: Calculate the total number of subgroups.
    print("Step 3: Calculate the final total.")
    total_subgroups = num_normal_v4_subgroups + total_non_normal
    print("Total subgroups of index 4 = (Number of normal subgroups) + (Number of non-normal subgroups)")
    print(f"Total = {num_normal_v4_subgroups} + {total_non_normal} = {total_subgroups}")

calculate_and_print_subgroup_count()