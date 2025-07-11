import sys

def analyze_crystal_classes():
    """
    Analyzes the 32 point crystal classes to find which are achiral,
    non-polar, and optically active.
    """
    # Definition of the 32 crystal classes with their properties.
    # A class is optically active if and only if it is chiral.
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

    print("Searching for crystal classes with the following properties:")
    print("1. Achiral (not chiral)")
    print("2. Non-polar")
    print("3. Optically Active (must be chiral)")
    print("-" * 30)

    found_classes = []
    for c in crystal_classes:
        # The condition for optical activity is chirality itself.
        is_optically_active = c['chiral']
        is_achiral = not c['chiral']
        is_non_polar = not c['polar']

        if is_achiral and is_non_polar and is_optically_active:
            found_classes.append(c['name'])

    if not found_classes:
        print("Result: No crystal classes were found.")
        print("\nExplanation:")
        print("The conditions are mutually exclusive. A crystal class is optically active if and only if it is chiral.")
        print("An achiral class, by definition, cannot be chiral, and therefore cannot be optically active.")
        print("Therefore, it is impossible for a crystal class to be both achiral and optically active.")
    else:
        # This part of the code is logically unreachable.
        print(f"Found matching crystal classes: {', '.join(found_classes)}")
    
    # Final answer as per instruction format.
    # Redirecting to stderr to not interfere with stdout for a potential programmatic user.
    print("\n<<<There are no crystal classes that fit the description.>>>", file=sys.stderr)


if __name__ == '__main__':
    analyze_crystal_classes()