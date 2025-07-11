def find_crystal_classes():
    """
    This script identifies crystal classes that are simultaneously achiral,
    non-polar, and optically active based on crystallographic rules.
    """

    # 1. Define classes based on their properties.

    # Optically active classes are those with a non-zero gyration tensor.
    # This includes the 11 chiral classes plus 6 specific achiral classes.
    optically_active_classes = {
        '1', '2', '222', '3', '32', '4', '422', '6', '622', '23', '432',  # Chiral
        'm', 'mm2', '-4', '3m', '4mm', '6mm'                               # Achiral
    }

    # Chiral classes lack any improper rotation axis (m, -1, -3, -4, -6).
    chiral_classes = {
        '1', '2', '222', '3', '32', '4', '422', '6', '622', '23', '432'
    }

    # Polar classes have a unique direction not related to its opposite by symmetry.
    polar_classes = {
        '1', '2', '3', '4', '6',
        'm', 'mm2', '3m', '4mm', '6mm'
    }

    # 2. Filter the classes step-by-step.

    print("Step 1: Identify optically active classes that are also achiral.")
    # An achiral class is any class not in the chiral_classes set.
    active_and_achiral = {c for c in optically_active_classes if c not in chiral_classes}
    print(f"Resulting classes: {sorted(list(active_and_achiral))}")

    print("\nStep 2: From the above set, identify the classes that are non-polar.")
    # A non-polar class is any class not in the polar_classes set.
    final_result = {c for c in active_and_achiral if c not in polar_classes}
    
    print("\nFinal Result:")
    print("The crystal class(es) that are achiral, non-polar, and optically active is/are:")
    if not final_result:
        print("None")
    else:
        # The prompt requires outputting each number/name in the final equation.
        # We will print each resulting class name.
        for crystal_class in sorted(list(final_result)):
            print(crystal_class)

find_crystal_classes()