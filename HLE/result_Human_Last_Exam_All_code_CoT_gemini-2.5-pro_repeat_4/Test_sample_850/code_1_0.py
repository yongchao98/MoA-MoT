def find_crystal_classes():
    """
    Analyzes the 32 crystal classes to find those that meet a set of
    contradictory criteria: achiral, non-polar, and optically active.
    """
    # Data for the 32 crystal classes (point groups).
    # A crystal is:
    # - Chiral if it lacks any roto-inversion axes (like mirror planes 'm' or
    #   a center of inversion 'i'). Chirality is the necessary condition for
    #   intrinsic optical activity.
    # - Polar if it has a unique direction whose two ends are not related by
    #   any symmetry operation in the point group.
    crystal_classes = [
        # Triclinic
        {'name': '1', 'is_chiral': True, 'is_polar': True},
        {'name': '-1', 'is_chiral': False, 'is_polar': False},
        # Monoclinic
        {'name': '2', 'is_chiral': True, 'is_polar': True},
        {'name': 'm', 'is_chiral': False, 'is_polar': True},
        {'name': '2/m', 'is_chiral': False, 'is_polar': False},
        # Orthorhombic
        {'name': '222', 'is_chiral': True, 'is_polar': False},
        {'name': 'mm2', 'is_chiral': False, 'is_polar': True},
        {'name': 'mmm', 'is_chiral': False, 'is_polar': False},
        # Tetragonal
        {'name': '4', 'is_chiral': True, 'is_polar': True},
        {'name': '-4', 'is_chiral': False, 'is_polar': False},
        {'name': '4/m', 'is_chiral': False, 'is_polar': False},
        {'name': '422', 'is_chiral': True, 'is_polar': False},
        {'name': '4mm', 'is_chiral': False, 'is_polar': True},
        {'name': '-42m', 'is_chiral': False, 'is_polar': False},
        {'name': '4/mmm', 'is_chiral': False, 'is_polar': False},
        # Trigonal
        {'name': '3', 'is_chiral': True, 'is_polar': True},
        {'name': '-3', 'is_chiral': False, 'is_polar': False},
        {'name': '32', 'is_chiral': True, 'is_polar': False},
        {'name': '3m', 'is_chiral': False, 'is_polar': True},
        {'name': '-3m', 'is_chiral': False, 'is_polar': False},
        # Hexagonal
        {'name': '6', 'is_chiral': True, 'is_polar': True},
        {'name': '-6', 'is_chiral': False, 'is_polar': False},
        {'name': '6/m', 'is_chiral': False, 'is_polar': False},
        {'name': '622', 'is_chiral': True, 'is_polar': False},
        {'name': '6mm', 'is_chiral': False, 'is_polar': True},
        {'name': '-6m2', 'is_chiral': False, 'is_polar': False},
        {'name': '6/mmm', 'is_chiral': False, 'is_polar': False},
        # Cubic
        {'name': '23', 'is_chiral': True, 'is_polar': False},
        {'name': 'm-3', 'is_chiral': False, 'is_polar': False},
        {'name': '432', 'is_chiral': True, 'is_polar': False},
        {'name': '-43m', 'is_chiral': False, 'is_polar': False},
        {'name': 'm-3m', 'is_chiral': False, 'is_polar': False},
    ]

    # The user requests classes that are:
    # 1. Achiral (is_chiral = False)
    # 2. Non-polar (is_polar = False)
    # 3. Optically active (which requires is_chiral = True)
    # These conditions are contradictory. The code will find an empty set.

    results = [
        c for c in crystal_classes
        if (not c['is_chiral']) and   # Achiral
           (not c['is_polar']) and    # Non-polar
           (c['is_chiral'])           # Optically Active (i.e., Chiral)
    ]

    # --- Output ---
    print("Searching for crystal classes with the following properties:")
    print("1. Achiral")
    print("2. Non-polar")
    print("3. Optically Active (which requires chirality)\n")
    print("As explained, conditions 1 and 3 are mutually exclusive.")
    print("----------------------------------------------------------")

    if not results:
        print("Result: No crystal classes were found that are simultaneously achiral, non-polar, and optically active.")
    else:
        # This code block is unreachable due to the logical contradiction
        print("Found matching crystal classes:")
        for c in results:
            print(f"- {c['name']}")

if __name__ == '__main__':
    find_crystal_classes()