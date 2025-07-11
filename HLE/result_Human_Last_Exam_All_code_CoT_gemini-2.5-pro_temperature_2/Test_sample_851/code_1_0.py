def find_optical_activity_classes():
    """
    Identifies crystal classes that are achiral, non-polar,
    and have symmetry allowing for directional optical activity (non-centrosymmetric).
    """
    # Properties for the 32 crystal classes (point groups).
    # Format: { 'Hermann-Mauguin symbol': [is_chiral, is_polar, is_centrosymmetric] }
    crystal_classes = {
        # Triclinic
        "1":        [True,  True,  False],
        "-1":       [False, False, True ],
        # Monoclinic
        "2":        [True,  True,  False],
        "m":        [False, True,  False],
        "2/m":      [False, False, True ],
        # Orthorhombic
        "222":      [True,  False, False],
        "mm2":      [False, True,  False],
        "mmm":      [False, False, True ],
        # Tetragonal
        "4":        [True,  True,  False],
        "-4":       [False, False, False],
        "4/m":      [False, False, True ],
        "422":      [True,  False, False],
        "4mm":      [False, True,  False],
        "-42m":     [False, False, False],
        "4/mmm":    [False, False, True ],
        # Trigonal
        "3":        [True,  True,  False],
        "-3":       [False, False, True ],
        "32":       [True,  False, False],
        "3m":       [False, True,  False],
        "-3m":      [False, False, True ],
        # Hexagonal
        "6":        [True,  True,  False],
        "-6":       [False, False, False],
        "6/m":      [False, False, True ],
        "622":      [True,  False, False],
        "6mm":      [False, True,  False],
        "-62m":     [False, False, False],
        "6/mmm":    [False, False, True ],
        # Cubic
        "23":       [True,  False, False],
        "m-3":      [False, False, True ],
        "432":      [True,  False, False],
        "-43m":     [False, False, False],
        "m-3m":     [False, False, True ],
    }

    print("Analyzing 32 crystal classes based on three conditions:")
    print("1. Achiral (not chiral)")
    print("2. Non-polar (not polar)")
    print("3. Allows directional optical activity (non-centrosymmetric)\n")

    candidate_classes = []
    for name, properties in crystal_classes.items():
        is_chiral, is_polar, is_centrosymmetric = properties
        
        # Apply the conditions:
        # Achiral AND Non-polar AND Non-centrosymmetric
        if not is_chiral and not is_polar and not is_centrosymmetric:
            candidate_classes.append(name)

    print("The crystal classes satisfying all conditions are:")
    # Print each identified class as requested
    for class_name in sorted(candidate_classes):
        print(class_name)

    print("\nComparing this list to the answer choices:")
    print("A. m and mm2 -> Incorrect (polar)")
    print("B. -6, -62m, and -43m -> Correct (subset of the full list)")
    print("C. 3m, 4m, and 6mm -> Incorrect (polar, and 4m is a typo for 4mm)")
    print("D. -4 and -42m -> Correct (subset of the full list)")
    print("E. 1, 2, 3, 4, and 6 -> Incorrect (chiral and polar)")
    print("\nBoth B and D are valid subsets. However, only one option can be chosen. Both are correct statements of fact. Option B is chosen as it contains a more diverse set of crystal systems (hexagonal and cubic).")

if __name__ == '__main__':
    find_optical_activity_classes()