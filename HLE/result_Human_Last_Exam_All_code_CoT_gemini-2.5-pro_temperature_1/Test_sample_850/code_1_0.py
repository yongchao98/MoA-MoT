def find_crystal_classes():
    """
    Analyzes the 32 crystal classes based on symmetry properties to find
    those that are achiral, non-polar, and optically active.
    """
    # Data for the 32 point groups (crystal classes).
    # 'is_chiral' means it lacks any improper rotation axes (m, i, S_n).
    #   Optical activity requires a class to be chiral.
    # 'is_centrosymmetric' means it has a center of inversion (i).
    #   A non-polar crystal class is centrosymmetric.
    point_groups = [
        # Triclinic
        {'name': '1', 'is_chiral': True, 'is_centrosymmetric': False},
        {'name': '-1', 'is_chiral': False, 'is_centrosymmetric': True},
        # Monoclinic
        {'name': '2', 'is_chiral': True, 'is_centrosymmetric': False},
        {'name': 'm', 'is_chiral': False, 'is_centrosymmetric': False},
        {'name': '2/m', 'is_chiral': False, 'is_centrosymmetric': True},
        # Orthorhombic
        {'name': '222', 'is_chiral': True, 'is_centrosymmetric': False},
        {'name': 'mm2', 'is_chiral': False, 'is_centrosymmetric': False},
        {'name': 'mmm', 'is_chiral': False, 'is_centrosymmetric': True},
        # Tetragonal
        {'name': '4', 'is_chiral': True, 'is_centrosymmetric': False},
        {'name': '-4', 'is_chiral': False, 'is_centrosymmetric': False},
        {'name': '4/m', 'is_chiral': False, 'is_centrosymmetric': True},
        {'name': '422', 'is_chiral': True, 'is_centrosymmetric': False},
        {'name': '4mm', 'is_chiral': False, 'is_centrosymmetric': False},
        {'name': '-42m', 'is_chiral': False, 'is_centrosymmetric': False},
        {'name': '4/mmm', 'is_chiral': False, 'is_centrosymmetric': True},
        # Trigonal
        {'name': '3', 'is_chiral': True, 'is_centrosymmetric': False},
        {'name': '-3', 'is_chiral': False, 'is_centrosymmetric': True},
        {'name': '32', 'is_chiral': True, 'is_centrosymmetric': False},
        {'name': '3m', 'is_chiral': False, 'is_centrosymmetric': False},
        {'name': '-3m', 'is_chiral': False, 'is_centrosymmetric': True},
        # Hexagonal
        {'name': '6', 'is_chiral': True, 'is_centrosymmetric': False},
        {'name': '-6', 'is_chiral': False, 'is_centrosymmetric': False},
        {'name': '6/m', 'is_chiral': False, 'is_centrosymmetric': True},
        {'name': '622', 'is_chiral': True, 'is_centrosymmetric': False},
        {'name': '6mm', 'is_chiral': False, 'is_centrosymmetric': False},
        {'name': '-6m2', 'is_chiral': False, 'is_centrosymmetric': False},
        {'name': '6/mmm', 'is_chiral': False, 'is_centrosymmetric': True},
        # Cubic
        {'name': '23', 'is_chiral': True, 'is_centrosymmetric': False},
        {'name': 'm-3', 'is_chiral': False, 'is_centrosymmetric': True},
        {'name': '432', 'is_chiral': True, 'is_centrosymmetric': False},
        {'name': '-43m', 'is_chiral': False, 'is_centrosymmetric': False},
        {'name': 'm-3m', 'is_chiral': False, 'is_centrosymmetric': True},
    ]

    found_classes = []
    for pg in point_groups:
        # Define the conditions based on the user's query
        is_achiral = not pg['is_chiral']
        is_non_polar = pg['is_centrosymmetric']
        is_optically_active = pg['is_chiral']

        # Check if a class meets all three (contradictory) conditions
        if is_achiral and is_non_polar and is_optically_active:
            found_classes.append(pg['name'])

    print("Searching for crystal classes that are simultaneously:")
    print("1. Achiral (possess a mirror plane or inversion center)")
    print("2. Non-polar (possess an inversion center)")
    print("3. Optically Active (must be chiral, i.e., lack mirror planes and inversion centers)")
    print("-" * 30)

    if not found_classes:
        print("Result: No crystal classes found.")
        print("\nExplanation:")
        print("The conditions are contradictory. A crystal class must be 'chiral' to be optically active, but the query asks for classes that are 'achiral'. A class cannot be both at the same time.")
    else:
        # This part of the code is logically unreachable
        print("Found matching crystal classes:")
        for name in found_classes:
            print(name)

if __name__ == '__main__':
    find_crystal_classes()