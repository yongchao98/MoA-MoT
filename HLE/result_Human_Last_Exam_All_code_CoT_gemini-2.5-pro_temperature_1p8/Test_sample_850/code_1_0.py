def solve_crystal_class_query():
    """
    Identifies crystal classes that are achiral, non-polar, and optically active,
    and explains the relationship between these properties.
    """
    # The 32 crystallographic point groups (crystal classes)
    all_32_classes = {
        '1', '-1', '2', 'm', '2/m', '222', 'mm2', 'mmm', '4', '-4', '4/m',
        '422', '4mm', '-42m', '4/mmm', '3', '-3', '32', '3m', '-3m', '6',
        '-6', '6/m', '622', '6mm', '-62m', '6/mmm', '23', 'm-3', '432',
        '-43m', 'm-3m'
    }

    # CRITERION 1: Optically Active
    # Optical activity is only possible in chiral crystal classes, which lack any
    # roto-inversion axis (including mirror planes 'm' and inversion centers '-1').
    optically_active_classes = {
        '1', '2', '222', '3', '32', '4', '422', '6', '622', '23', '432'
    }

    # CRITERION 2: Achiral
    # Achiral classes are the opposite of chiral classes. They are the set of all
    # classes that are NOT optically active.
    achiral_classes = all_32_classes - optically_active_classes

    # CRITERION 3: Non-polar
    # Polar classes exhibit a spontaneous electric dipole. Non-polar classes do not.
    # There are 10 polar classes. All others are non-polar.
    polar_classes = {'1', '2', 'm', 'mm2', '3', '3m', '4', '4mm', '6', '6mm'}
    non_polar_classes = all_32_classes - polar_classes

    # Find the intersection of all three sets of criteria.
    # Equation: Result = (Achiral Classes) ∩ (Non-polar Classes) ∩ (Optically Active Classes)
    result_set = achiral_classes.intersection(non_polar_classes).intersection(optically_active_classes)

    print("--- Analysis of Crystal Class Properties ---")
    print("The query asks for crystal classes that are simultaneously:")
    print("  1. Achiral")
    print("  2. Non-polar")
    print("  3. Optically Active\n")

    print("There is a fundamental conflict in these requirements:")
    print("  - Optical Activity REQUIRES a class to be CHIRAL.")
    print("  - A class cannot be both CHIRAL and ACHIRAL at the same time.\n")

    print("Therefore, no crystal class can be both Achiral and Optically Active.")
    print("The intersection of these sets must be empty.\n")
    print("--- Calculation ---")
    print(f"Set of Achiral Classes ({len(achiral_classes)} members):")
    print(f"{sorted(list(achiral_classes))}\n")
    print(f"Set of Non-polar Classes ({len(non_polar_classes)} members):")
    print(f"{sorted(list(non_polar_classes))}\n")
    print(f"Set of Optically Active Classes ({len(optically_active_classes)} members):")
    print(f"{sorted(list(optically_active_classes))}\n")

    print("--- Final Equation and Result ---")
    print("Result = (Achiral ∩ Non-polar ∩ Optically Active)")
    # The "final equation" prints the members of the final result set.
    # In this case, there are none.
    print(f"Number of matching crystal classes: {len(result_set)}")
    if not result_set:
        print("Final Result Set: {} (The empty set)")
    else:
        print(f"Final Result Set: {sorted(list(result_set))}")

solve_crystal_class_query()