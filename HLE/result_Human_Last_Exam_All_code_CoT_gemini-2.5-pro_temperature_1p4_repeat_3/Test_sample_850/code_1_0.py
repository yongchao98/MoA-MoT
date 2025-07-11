def find_optically_active_achiral_classes():
    """
    Analyzes the 32 crystal classes to find which are achiral, non-polar,
    and have the correct symmetry for optical activity.
    """

    # The 11 chiral (and thus optically active) crystal classes.
    # Optical activity is possible ONLY in these classes.
    optically_active_classes = {
        '1', '2', '3', '4', '6', '222', '32', '422', '622', '23', '432'
    }

    # The 32 crystallographic point groups (crystal classes).
    all_crystal_classes = {
        '1', '-1', '2', 'm', '2/m', '222', 'mm2', 'mmm', '4', '-4', '4/m',
        '422', '4mm', '-42m', '4/mmm', '3', '-3', '32', '3m', '-3m', '6',
        '-6', '6/m', '622', '6mm', '-62m', '6/mmm', '23', 'm-3', '432',
        '-43m', 'm-3m'
    }

    # Achiral classes are all classes that are NOT chiral (optically active).
    # They contain a center of inversion or a mirror plane.
    achiral_classes = all_crystal_classes - optically_active_classes

    # The question asks for crystal classes that are BOTH achiral AND optically active.
    # We find the intersection of these two sets.
    result_set = achiral_classes.intersection(optically_active_classes)

    print("Step 1: Define the condition for optical activity.")
    print("A crystal class must be CHIRAL to be optically active.")
    print("The 11 chiral (optically active) classes are:")
    print(sorted(list(optically_active_classes)))
    print("-" * 20)

    print("Step 2: Define the condition for being achiral.")
    print("An achiral class is any class that is NOT chiral.")
    print("The 21 achiral classes are:")
    print(sorted(list(achiral_classes)))
    print("-" * 20)

    print("Step 3: Find the crystal classes that satisfy BOTH conditions.")
    print("We are looking for classes in the intersection of the 'achiral' set and the 'optically active' set.")
    print(f"Intersection of Achiral and Optically Active sets: {list(result_set)}")
    print("-" * 20)

    print("Conclusion:")
    print("There is a fundamental contradiction in the request. A crystal cannot be both achiral and optically active.")
    print("As shown by the empty intersection set, no such crystal classes exist.")

find_optically_active_achiral_classes()