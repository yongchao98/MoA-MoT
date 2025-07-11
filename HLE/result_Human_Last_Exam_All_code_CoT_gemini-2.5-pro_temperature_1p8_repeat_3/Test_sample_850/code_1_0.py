def solve_crystal_class_query():
    """
    Analyzes crystallographic point groups to find those matching specific symmetry criteria.
    """

    # The 11 chiral (and thus optically active) crystal classes.
    # These lack a center of inversion or a mirror plane.
    chiral_classes = {'1', '2', '222', '3', '32', '4', '422', '6', '622', '23', '432'}

    # The 21 achiral crystal classes are all other classes.
    # We can define them for completeness.
    achiral_classes = {
        '-1', 'm', '2/m',
        'mmm', 'mm2',
        '-3', '3m', '-3m',
        '4/m', '-4', '4mm', '4/mmm', '42m',
        '6/m', '-6', '6mm', '6/mmm', '6m2',
        'm-3', '43m', 'm-3m'
    }

    # The 10 polar crystal classes.
    # These have a unique polar axis, leading to properties like pyroelectricity.
    polar_classes = {'1', '2', '3', '4', '6', 'm', 'mm2', '3m', '4mm', '6mm'}

    print("Analyzing the user's original query: 'achiral', 'non-polar', and 'optically active'")
    print("-" * 70)
    # A class is optically active if and only if it is chiral.
    optically_active_classes = chiral_classes

    # The query asks for classes that are in the intersection of 'achiral_classes' and 'optically_active_classes'.
    contradiction_set = achiral_classes.intersection(optically_active_classes)

    print(f"A crystal class is optically active if and only if it is chiral.")
    print(f"By definition, a class cannot be both chiral and achiral.")
    print(f"The set of classes that are both 'achiral' and 'optically active' is: {contradiction_set if contradiction_set else 'Empty'}\n")


    print("Assuming the intended query was for 'chiral' (optically active) and 'non-polar' classes:")
    print("-" * 70)
    # We find the set of chiral classes that are not polar.
    chiral_non_polar_classes = chiral_classes.difference(polar_classes)

    # Sort the result for clear presentation
    sorted_result = sorted(list(chiral_non_polar_classes), key=lambda x: (len(x), x))

    print("The optically active (chiral) and non-polar crystal classes are:")
    # Print the final list
    print(', '.join(sorted_result))


solve_crystal_class_query()
<<<222, 32, 422, 622, 23, 432>>>