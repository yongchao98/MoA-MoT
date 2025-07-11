def solve_isomer_count():
    """
    Calculates the number of isomers formed from the reaction of
    2,5-di(2-pyridyl)thiazolo[5,4-d]thiazole with
    cis-dichlorobis(bipyridine)ruthenium(II).

    The product is assumed to be [Ru(bpy)(py-tz)(Cl)2], which has the
    general form [M(AA)(BC)X2].
    """
    print("Step 1: Counting isomers for the trans-dichloro configuration.")
    # In the trans-[M(AA)(BC)X2] complex, two different orientations of the
    # unsymmetrical BC ligand relative to the AA ligand create two diastereomers.
    # Each diastereomer is chiral, existing as an enantiomeric pair.
    num_trans_diastereomers = 2
    enantiomers_per_pair = 2
    num_trans_isomers = num_trans_diastereomers * enantiomers_per_pair
    print(f"For the trans-Cl2 case, there are {num_trans_diastereomers} diastereomers, each with an enantiomer.")
    print(f"Number of trans isomers = {num_trans_diastereomers} * {enantiomers_per_pair} = {num_trans_isomers}")
    print("-" * 20)

    print("Step 2: Counting isomers for the cis-dichloro configuration.")

    # Subcase 2a: 'fac' arrangement of chelates (bidentate ligands are cis)
    # This arrangement forms a single chiral structure, i.e., one pair of enantiomers.
    num_fac_pairs = 1
    num_fac_isomers = num_fac_pairs * enantiomers_per_pair
    print("For the cis-Cl2, fac-chelate case, there is 1 pair of enantiomers.")
    print(f"Number of cis-fac isomers = {num_fac_pairs} * {enantiomers_per_pair} = {num_fac_isomers}")

    # Subcase 2b: 'mer' arrangement of chelates (bidentate ligands are trans)
    # This arrangement forms two distinct chiral diastereomers.
    num_mer_diastereomers = 2
    num_mer_isomers = num_mer_diastereomers * enantiomers_per_pair
    print("For the cis-Cl2, mer-chelate case, there are 2 diastereomers, each with an enantiomer.")
    print(f"Number of cis-mer isomers = {num_mer_diastereomers} * {enantiomers_per_pair} = {num_mer_isomers}")
    print("-" * 20)

    print("Step 3: Calculating the total number of isomers.")
    total_isomers = num_trans_isomers + num_fac_isomers + num_mer_isomers
    print(f"Total isomers = (trans isomers) + (cis-fac isomers) + (cis-mer isomers)")
    print(f"Total isomers = {num_trans_isomers} + {num_fac_isomers} + {num_mer_isomers}")
    print(f"Total = {total_isomers}")

solve_isomer_count()