def calculate_mercedesbenzene_carbons():
    """
    Calculates the number of carbon atoms in the fictitious molecule
    "mercedesbenzene", interpreted as the chemical compound triphenylene.
    """
    # A central benzene ring has 6 carbon atoms.
    central_ring = 6
    
    # The three "points" of the Mercedes star are interpreted as three
    # additional benzene rings fused to the central ring.
    # When a ring is fused, it shares 2 carbons, so it adds 4 new carbons.
    first_fused_ring_added_carbons = 4
    second_fused_ring_added_carbons = 4
    third_fused_ring_added_carbons = 4

    # Sum the carbons to find the total.
    total_carbons = central_ring + first_fused_ring_added_carbons + second_fused_ring_added_carbons + third_fused_ring_added_carbons

    # Print the final equation showing each number.
    print("The molecule is interpreted as triphenylene, which consists of a central ring and three fused rings.")
    print("The total number of carbons is calculated as follows:")
    print(f"{central_ring} + {first_fused_ring_added_carbons} + {second_fused_ring_added_carbons} + {third_fused_ring_added_carbons} = {total_carbons}")

calculate_mercedesbenzene_carbons()