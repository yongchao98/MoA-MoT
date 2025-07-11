def find_crystal_classes():
    """
    Defines all 32 crystal classes and their properties, then searches
    for classes matching the user's contradictory criteria.
    """
    # Data for the 32 crystallographic point groups (crystal classes)
    # 'is_chiral': Lacks inversion center and mirror planes.
    # 'is_polar': Has a unique directional axis.
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

    # The fundamental condition for optical activity is chirality.
    # A crystal class is optically active if and only if it is chiral.
    def is_optically_active(crystal_class):
        return crystal_class['is_chiral']

    # Find classes that are achiral, non-polar, AND optically active
    results = []
    for cc in crystal_classes:
        # The user's query conditions:
        is_achiral_req = not cc['is_chiral']
        is_non_polar_req = not cc['is_polar']
        is_optically_active_req = is_optically_active(cc)

        if is_achiral_req and is_non_polar_req and is_optically_active_req:
            results.append(cc['name'])

    # Print the results and explanation
    print("Searching for crystal classes that are simultaneously:")
    print("1. Achiral (not chiral)")
    print("2. Non-polar")
    print("3. Optically active")
    print("-" * 30)

    if not results:
        print("Result: No such crystal classes exist.")
        print("\nExplanation:")
        print("The requirement for a crystal to be optically active is that it must be chiral.")
        print("The query asks for classes that are 'achiral' and 'optically active', which is a contradiction in terms.")
        print("Therefore, the set of crystal classes that meet these criteria is empty.")
    else:
        # This part of the code will not be reached, but is included for completeness.
        print("Found matching crystal classes:")
        for name in results:
            print(name)

if __name__ == '__main__':
    find_crystal_classes()