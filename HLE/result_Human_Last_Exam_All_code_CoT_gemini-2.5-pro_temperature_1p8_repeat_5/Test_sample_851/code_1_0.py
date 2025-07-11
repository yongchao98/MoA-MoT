def solve_crystal_class_problem():
    """
    Identifies crystal classes that are achiral, non-polar, and can exhibit optical activity.
    """
    # Dictionary of the 32 crystallographic point groups and their properties.
    # Format: { 'Hermann-Mauguin symbol': [is_centrosymmetric, is_achiral, is_polar] }
    # is_achiral is True if it's NOT in the 11 chiral (enantiomorphic) classes.
    point_groups = {
        # Triclinic
        "1":    [False, False, True],
        "-1":   [True,  True,  False],
        # Monoclinic
        "2":    [False, False, True],
        "m":    [False, True,  True],
        "2/m":  [True,  True,  False],
        # Orthorhombic
        "222":  [False, False, False],
        "mm2":  [False, True,  True],
        "mmm":  [True,  True,  False],
        # Tetragonal
        "4":    [False, False, True],
        "-4":   [False, True,  False],
        "4/m":  [True,  True,  False],
        "422":  [False, False, False],
        "4mm":  [False, True,  True],
        "-42m": [False, True,  False],
        "4/mmm":[True,  True,  False],
        # Trigonal
        "3":    [False, False, True],
        "-3":   [True,  True,  False],
        "32":   [False, False, False],
        "3m":   [False, True,  True],
        "-3m":  [True,  True,  False],
        # Hexagonal
        "6":    [False, False, True],
        "-6":   [False, True,  False],
        "6/m":  [True,  True,  False],
        "622":  [False, False, False],
        "6mm":  [False, True,  True],
        "-62m": [False, True,  False],
        "6/mmm":[True,  True,  False],
        # Cubic
        "23":   [False, False, False],
        "m-3":  [True,  True,  False],
        "432":  [False, False, False],
        "-43m": [False, True,  False],
        "m-3m": [True,  True,  False]
    }

    # Criteria for filtering:
    # 1. Can be optically active -> non-centrosymmetric
    # 2. Achiral
    # 3. Non-polar
    
    print("Applying the following criteria to the 32 point groups:")
    print("1. Achiral: The point group must possess an improper rotation axis (e.g., 'm' or '-n').")
    print("2. Non-polar: The point group must not be one of the 10 polar classes.")
    print("3. Capable of Optical Activity: The point group must be non-centrosymmetric.")
    print("-" * 30)

    result_classes = []
    for pg, props in point_groups.items():
        is_centrosymmetric, is_achiral, is_polar = props
        if not is_centrosymmetric and is_achiral and not is_polar:
            result_classes.append(pg)

    print("The crystal classes satisfying all three conditions are:")
    print(sorted(result_classes))
    print("-" * 30)

    print("Analyzing the given answer choices based on this result:")
    print("A. m and mm2: Incorrect. These classes are polar.")
    print("B. -6, -62m, and -43m: Correct. All three classes are in the list of results.")
    print("C. 3m, 4mm, and 6mm: Incorrect. These classes are polar.")
    print("D. -4 and -42m: Correct. Both classes are in the list of results.")
    print("E. 1, 2, 3, 4, and 6: Incorrect. These classes are chiral (not achiral) and polar.")
    print("\nThe complete set of classes is {'-4', '-42m', '-6', '-62m', '-43m'}. Both B and D are valid but incomplete lists. However, a choice must be made. Choice B contains more classes from multiple crystal systems.")

solve_crystal_class_problem()