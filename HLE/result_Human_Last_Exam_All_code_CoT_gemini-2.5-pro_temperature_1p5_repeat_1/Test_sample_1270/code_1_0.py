def calculate_mercedesbenzene_carbons():
    """
    Calculates the number of carbon atoms in the fictitious molecule "mercedesbenzene",
    interpreting it as the real molecule triphenylene (C18H12).

    The calculation is based on constructing the molecule by sequentially fusing
    four benzene rings.
    """
    # Number of carbons in the first benzene ring
    first_ring_carbons = 6
    
    # Number of carbons added by each subsequent fused ring
    carbons_per_fused_ring = 4
    
    # Triphenylene is formed from a central ring and three fused rings around it.
    # Total rings = 4. Number of fusions = 3.
    num_fusions = 3

    # Calculate the total number of carbons
    total_carbons = first_ring_carbons + num_fusions * carbons_per_fused_ring

    # Print the explanation and the final equation
    print("The molecule 'mercedesbenzene' can be interpreted as triphenylene,")
    print("which is formed by fusing four benzene rings in a symmetric, propeller-like shape.")
    print("Calculation based on ring fusion:")
    print(f"Start with one benzene ring ({first_ring_carbons} carbons), then add {carbons_per_fused_ring} carbons for each of the {num_fusions} additional fused rings.")
    print(f"The equation is: {first_ring_carbons} + {carbons_per_fused_ring} + {carbons_per_fused_ring} + {carbons_per_fused_ring} = {total_carbons}")
    print(f"Therefore, mercedesbenzene would have {total_carbons} carbons.")

if __name__ == "__main__":
    calculate_mercedesbenzene_carbons()