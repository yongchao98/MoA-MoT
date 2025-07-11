def solve_isomer_problem():
    """
    Calculates the number of isomers formed in the reaction of
    2,5-di(2-pyridyl)thiazolo[5,4-d]thiazole with cis-dichlorobis(bipyridine)ruthenium(II).
    """

    # Step 1: Analyze the product complex.
    # The reaction is a substitution of two chloride ligands (Cl-) by the bidentate ligand L
    # (L = 2,5-di(2-pyridyl)thiazolo[5,4-d]thiazole).
    # The starting complex is cis-[Ru(bpy)2Cl2].
    # The product complex is [Ru(bpy)2(L)]^2+.
    # This product is an octahedral tris-chelate complex of the form [M(A-A)2(B-B)],
    # where A-A is bipyridine (bpy) and B-B is the new ligand (L).
    # Both bpy and L are symmetric bidentate ligands.

    # Step 2: Determine the number of geometric isomers.
    # For a complex of the type [M(A-A)2(B-B)] with symmetric bidentate ligands,
    # there are two possible spatial arrangements for the ligands. These are diastereomers.
    # One has C2 symmetry, and the other has C1 symmetry.
    num_geometric_isomers = 2
    print(f"Number of geometric isomers (diastereomers): {num_geometric_isomers}")

    # Step 3: Determine the number of optical isomers (enantiomers) for each geometric isomer.
    # All octahedral complexes with three bidentate ligands are chiral (lack a mirror plane).
    # This means each geometric isomer exists as a pair of non-superimposable mirror images (enantiomers).
    enantiomers_per_geometric_isomer = 2
    print(f"Number of optical isomers per geometric isomer: {enantiomers_per_geometric_isomer}")

    # Step 4: Calculate the total number of isomers.
    # The total number of stereoisomers is the product of the number of geometric isomers
    # and the number of optical isomers per geometric isomer.
    total_isomers = num_geometric_isomers * enantiomers_per_geometric_isomer

    print("\nFinal Calculation:")
    print(f"The total number of isomers is {num_geometric_isomers} (geometric isomers) * {enantiomers_per_geometric_isomer} (enantiomers per pair).")
    print(f"Total isomers = {total_isomers}")

solve_isomer_problem()
