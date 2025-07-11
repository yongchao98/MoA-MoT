def calculate_mercedesbenzene_carbons():
    """
    Calculates the number of carbon atoms in the fictitious molecule "mercedesbenzene".

    This structure is interpreted as a central benzene ring with three additional
    benzene rings fused to it symmetrically, resembling the Mercedes-Benz logo.
    This structure is chemically known as triphenylene.
    """
    # Number of carbons in the central benzene ring
    central_ring_carbons = 6

    # Number of benzene rings fused to the central ring
    num_attached_rings = 3

    # Each attached benzene ring (6 carbons) shares 2 carbons with the central ring,
    # so it adds (6 - 2) = 4 new carbons to the molecule.
    new_carbons_per_ring = 4

    # Calculate the total number of carbon atoms
    total_carbons = central_ring_carbons + num_attached_rings * new_carbons_per_ring

    print("This molecule is envisioned as a central benzene ring with three other rings fused to it.")
    print(f"Carbons in the central ring: {central_ring_carbons}")
    print(f"Number of attached rings: {num_attached_rings}")
    print(f"New carbons added by each attached ring: {new_carbons_per_ring}")
    print("\nFinal Equation:")
    print(f"Total Carbons = {central_ring_carbons} + {num_attached_rings} * {new_carbons_per_ring} = {total_carbons}")

calculate_mercedesbenzene_carbons()