def find_crystal_classes():
    """
    Analyzes the 32 crystallographic point groups to find those that are
    achiral, non-polar, and optically active, and explains the result.
    """
    # Using standard Hermann-Mauguin notation.
    # Note: -n is used to represent n-bar (e.g., -1 is 1-bar).
    all_32_classes = {
        '1', '-1', '2', 'm', '2/m', '222', 'mm2', 'mmm',
        '4', '-4', '4/m', '422', '4mm', '-42m', '4/mmm',
        '3', '-3', '32', '3m', '-3m',
        '6', '-6', '6/m', '622', '6mm', '-6m2', '6/mmm',
        '23', 'm-3', '432', '-43m', 'm-3m'
    }

    # Condition 1: Optically Active (the 11 Chiral classes)
    # These lack any mirror planes (m) or inversion centers (-1).
    optically_active_classes = {
        '1', '2', '222', '3', '32', '4', '422', '6', '622', '23', '432'
    }

    # Condition 2: Achiral classes
    # These are all classes that are NOT chiral.
    achiral_classes = all_32_classes - optically_active_classes

    # Condition 3: Non-polar classes
    # Polar classes are {1, 2, m, mm2, 3, 3m, 4, 4mm, 6, 6mm}. Non-polar are all others.
    polar_classes = {'1', '2', 'm', 'mm2', '3', '3m', '4', '4mm', '6', '6mm'}
    non_polar_classes = all_32_classes - polar_classes

    # Find the intersection of classes that are: Achiral AND Non-polar AND Optically Active
    result_set = achiral_classes.intersection(non_polar_classes).intersection(optically_active_classes)
    
    # --- Explanation ---
    print("--- Analysis of Crystal Class Properties ---")
    print("The query seeks crystal classes that are simultaneously achiral, non-polar, and optically active.")
    print("\nFundamental Principle: A crystal class can only be optically active if it is chiral.")
    print(" - Chiral means it lacks mirror planes and inversion centers.")
    print(" - Achiral means it possesses a mirror plane or an inversion center.")
    print("\nConclusion: The conditions 'achiral' and 'optically active' are mutually exclusive. Therefore, no such crystal classes can exist.")
    print("The script below proves this by finding the intersection of these defined sets.")
    print("-" * 40)

    # --- Calculation ---
    print(f"Set 1: Achiral Classes = {sorted(list(achiral_classes))}")
    print(f"Set 2: Non-Polar Classes = {sorted(list(non_polar_classes))}")
    print(f"Set 3: Optically Active Classes = {sorted(list(optically_active_classes))}")
    
    # --- Final Result ---
    print("\nIntersection of (Achiral AND Non-Polar AND Optically Active):")
    if not result_set:
        print("The resulting set is empty.")
        print("\nFinal Answer: There are 0 crystal classes that satisfy all conditions.")
    else:
        # This branch will not be hit, but is included for logical completeness
        print(f"Resulting Classes: {sorted(list(result_set))}")


if __name__ == '__main__':
    find_crystal_classes()
