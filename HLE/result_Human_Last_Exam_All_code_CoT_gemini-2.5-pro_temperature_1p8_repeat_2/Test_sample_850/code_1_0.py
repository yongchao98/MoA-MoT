def solve_crystal_class_query():
    """
    Identifies crystal classes based on symmetry properties related to optical activity.
    """
    # The 32 crystallographic point groups (crystal classes) in Hermann-Mauguin notation
    all_32_classes = {
        '1', '-1', '2', 'm', '2/m', '222', 'mm2', 'mmm', '3', '-3', '32',
        '3m', '-3m', '4', '-4', '4/m', '422', '4mm', '-42m', '4/mmm', '6',
        '-6', '6/m', '622', '6mm', '-6m2', '6/mmm', '23', 'm-3', '432',
        '-43m', 'm-3m'
    }

    # A crystal class must be chiral (lacking rotoinversion axes like mirror planes or inversion centers) to be optically active.
    # There are 11 such classes.
    optically_active_classes = {'1', '2', '3', '4', '6', '222', '32', '422', '622', '23', '432'}
    chiral_classes = optically_active_classes

    # Achiral classes are all classes that are not chiral.
    achiral_classes = all_32_classes - chiral_classes

    # Polar classes are those with a unique axis of polarity.
    polar_classes = {'1', '2', 'm', 'mm2', '3', '3m', '4', '4mm', '6', '6mm'}

    # Non-polar classes are all classes that are not polar.
    non_polar_classes = all_32_classes - polar_classes

    print("--- Analysis of Crystal Classes for Optical Activity ---")
    print("\n1. Evaluating the user's specific query:")
    print("Query: Find crystal classes that are ACHIRAL, NON-POLAR, and OPTICALLY ACTIVE.\n")
    
    print("A fundamental principle of crystallography is that optical activity can only occur in CHIRAL crystal classes.")
    print("An ACHIRAL class, by definition, cannot be optically active. Therefore, we expect no results for this query.")
    
    # Calculate the intersection based on the user's query
    result_set = achiral_classes.intersection(non_polar_classes, optically_active_classes)
    
    print("\nResult of the search for {Achiral, Non-Polar, Optically Active} classes:")
    if not result_set:
        print("None")
    else:
        # This part of the code should not be reached
        print(sorted(list(result_set)))

    print("\n-------------------------------------------------------------")
    print("\n2. Clarification for a related, valid query:")
    print("Query: Find crystal classes that are CHIRAL, NON-POLAR, and OPTICALLY ACTIVE.\n")
    print("This query identifies the non-polar classes that do exhibit optical activity.")
    
    # Calculate the intersection for the clarified query
    clarified_result_set = chiral_classes.intersection(non_polar_classes)

    print("\nOptically active (chiral) classes that are also non-polar:")
    # We sort the list for consistent and readable output
    print(sorted(list(clarified_result_set)))

solve_crystal_class_query()