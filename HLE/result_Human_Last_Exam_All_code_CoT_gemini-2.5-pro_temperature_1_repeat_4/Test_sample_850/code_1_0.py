def find_crystal_classes():
    """
    Analyzes the 32 crystallographic point groups to find those that are
    achiral, non-polar, and optically active.
    """
    # Data for the 32 crystal classes (point groups).
    # 'chiral' implies optical activity.
    # 'polar' refers to the presence of at least one unique polar direction.
    crystal_classes = [
        # Triclinic
        {'name': '1', 'chiral': True, 'polar': True},
        {'name': '-1', 'chiral': False, 'polar': False},
        # Monoclinic
        {'name': '2', 'chiral': True, 'polar': True},
        {'name': 'm', 'chiral': False, 'polar': True},
        {'name': '2/m', 'chiral': False, 'polar': False},
        # Orthorhombic
        {'name': '222', 'chiral': True, 'polar': False},
        {'name': 'mm2', 'chiral': False, 'polar': True},
        {'name': 'mmm', 'chiral': False, 'polar': False},
        # Tetragonal
        {'name': '4', 'chiral': True, 'polar': True},
        {'name': '-4', 'chiral': False, 'polar': False},
        {'name': '4/m', 'chiral': False, 'polar': False},
        {'name': '422', 'chiral': True, 'polar': False},
        {'name': '4mm', 'chiral': False, 'polar': True},
        {'name': '-42m', 'chiral': False, 'polar': False},
        {'name': '4/mmm', 'chiral': False, 'polar': False},
        # Trigonal
        {'name': '3', 'chiral': True, 'polar': True},
        {'name': '-3', 'chiral': False, 'polar': False},
        {'name': '32', 'chiral': True, 'polar': False},
        {'name': '3m', 'chiral': False, 'polar': True},
        {'name': '-3m', 'chiral': False, 'polar': False},
        # Hexagonal
        {'name': '6', 'chiral': True, 'polar': True},
        {'name': '-6', 'chiral': False, 'polar': False},
        {'name': '6/m', 'chiral': False, 'polar': False},
        {'name': '622', 'chiral': True, 'polar': False},
        {'name': '6mm', 'chiral': False, 'polar': True},
        {'name': '-6m2', 'chiral': False, 'polar': False},
        {'name': '6/mmm', 'chiral': False, 'polar': False},
        # Cubic
        {'name': '23', 'chiral': True, 'polar': False},
        {'name': 'm-3', 'chiral': False, 'polar': False},
        {'name': '432', 'chiral': True, 'polar': False},
        {'name': '-43m', 'chiral': False, 'polar': False},
        {'name': 'm-3m', 'chiral': False, 'polar': False},
    ]

    print("--- Analysis of Crystal Class Properties ---")
    print("A crystal class must be CHIRAL to be optically active.")
    print("A chiral class, by definition, lacks mirror planes and inversion centers.")
    print("An ACHIRAL class, by definition, possesses a mirror plane or an inversion center.")
    print("\nConclusion: The conditions 'achiral' and 'optically active' are mutually exclusive.")
    print("Therefore, no crystal class can satisfy both conditions.\n")

    # Search for classes matching the user's criteria
    found_classes = []
    for cc in crystal_classes:
        is_achiral = not cc['chiral']
        is_non_polar = not cc['polar']
        is_optically_active = cc['chiral'] # Optical activity requires chirality

        if is_achiral and is_non_polar and is_optically_active:
            found_classes.append(cc['name'])

    print("--- Search Results ---")
    print("Searching for crystal classes that are simultaneously:")
    print("1. Achiral")
    print("2. Non-polar")
    print("3. Optically active")
    print("\nMatching crystal classes found:")
    print(found_classes)

if __name__ == '__main__':
    find_crystal_classes()