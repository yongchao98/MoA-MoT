def explain_crystal_class_properties():
    """
    This function analyzes the user's request for crystal classes that are
    achiral, non-polar, and optically active, and explains the inherent
    contradiction in the query.
    """
    # A crystal is optically active if and only if its point group is chiral.
    # Chiral point groups lack mirror planes (m) and inversion centers (i).
    optically_active_classes = {
        '1', '2', '3', '4', '6', '222', '32', '422', '622', '23', '432'
    }

    # An achiral class is one that is NOT chiral. By definition, it must contain
    # a mirror plane or an inversion center.
    all_32_point_groups = {
        '1', '-1', '2', 'm', '2/m', '222', 'mm2', 'mmm',
        '3', '-3', '32', '3m', '-3m',
        '4', '-4', '4/m', '422', '4mm', '-42m', '4/mmm',
        '6', '-6', '6/m', '622', '6mm', '-62m', '6/mmm',
        '23', 'm-3', '432', '-43m', 'm-3m'
    }
    achiral_classes = all_32_point_groups - optically_active_classes

    # The request is for classes that are in BOTH the 'optically_active_classes' set
    # AND the 'achiral_classes' set. Let's find the intersection.
    intersection = optically_active_classes.intersection(achiral_classes)

    print("--- Analysis of the Request ---")
    print("\nThe request is for crystal classes that have the following properties:")
    print("1. Achiral")
    print("2. Non-polar")
    print("3. Optically Active")
    print("\n--- Examining the Core Conflict ---")
    print("\nCondition (3) 'Optically Active' requires a crystal's point group to be CHIRAL.")
    print("The 11 chiral (and thus optically active) point groups are:")
    print(f"-> {sorted(list(optically_active_classes))}")
    
    print("\nCondition (1) 'Achiral' requires a crystal's point group to be ACHIRAL (i.e., not chiral).")
    print("The 21 achiral point groups are:")
    print(f"-> {sorted(list(achiral_classes))}")
    
    print("\n--- Conclusion ---")
    print("A crystal class cannot be both CHIRAL and ACHIRAL at the same time.")
    print("The requirement for optical activity (chirality) and the condition of being achiral are mutually exclusive.")
    
    if not intersection:
        print("\nTherefore, there are no crystal classes that can satisfy these contradictory conditions.")
        print("The set of 'achiral crystal classes with symmetry for optical activity' is empty.")
    else:
        # This code is unreachable due to the definitions used
        print(f"\nError: A contradictory set was found: {intersection}")

explain_crystal_class_properties()