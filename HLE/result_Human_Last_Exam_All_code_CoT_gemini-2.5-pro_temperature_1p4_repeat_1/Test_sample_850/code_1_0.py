def find_crystal_classes():
    """
    Identifies crystal classes based on a combination of symmetry properties.

    The user is asking for crystal classes that are simultaneously:
    1. Achiral
    2. Non-polar
    3. Optically Active

    This script demonstrates that no such crystal classes exist because the conditions
    "achiral" and "optically active" are mutually exclusive.
    """

    # Optical activity is a property of chiral materials. Therefore, the set of
    # optically active classes is identical to the set of chiral classes.
    # These 11 classes lack any improper symmetry (inversion, mirror planes, rotoinversion).
    optically_active_classes = {
        '1', '2', '3', '4', '6', '222', '32', '422', '622', '23', '432'
    }

    # Achiral classes are, by definition, all classes that are NOT chiral.
    # These 21 classes all contain at least one improper symmetry element.
    all_classes = {
        '1', '-1', '2', 'm', '2/m', '222', 'mm2', 'mmm', '4', '-4', '4/m',
        '422', '4mm', '-42m', '4/mmm', '3', '-3', '32', '3m', '-3m', '6',
        '-6', '6/m', '622', '6mm', '-6m2', '6/mmm', '23', 'm-3', '432',
        '-43m', 'm-3m'
    }
    achiral_classes = all_classes - optically_active_classes

    # Non-polar classes lack a unique polar direction. They include all
    # centrosymmetric classes and some non-centrosymmetric ones.
    # The 10 polar classes are: 1, 2, 3, 4, 6, m, mm2, 3m, 4mm, 6mm
    polar_classes = {'1', '2', '3', '4', '6', 'm', 'mm2', '3m', '4mm', '6mm'}
    non_polar_classes = all_classes - polar_classes

    # Now, find the intersection of the three sets as requested by the user.
    # This corresponds to: (Achiral) AND (Non-polar) AND (Optically Active)
    result = achiral_classes.intersection(non_polar_classes).intersection(optically_active_classes)

    print("This script checks for crystal classes with three properties: Achiral, Non-polar, and Optically Active.\n")
    print("Principle: A crystal class is optically active if and only if it is chiral.")
    print("Principle: A crystal class is achiral if it is not chiral.")
    print("Conclusion: Therefore, 'optically active' and 'achiral' are mutually exclusive properties.\n")

    print("The sets of crystal classes based on each property are:")
    print(f"1. Optically Active Classes: {sorted(list(optically_active_classes))}")
    print(f"2. Achiral Classes: {sorted(list(achiral_classes))}")
    print(f"3. Non-polar Classes: {sorted(list(non_polar_classes))}\n")

    print("We now find the intersection of these three sets.")
    # The "equation" is the set intersection operation. We output each number (class name) in the result.
    print(f"(Achiral) \u2229 (Non-polar) \u2229 (Optically Active) = {list(result)}")
    print("\nThe resulting set is empty. This proves that there are no crystal classes that are simultaneously achiral, non-polar, and optically active.")

if __name__ == '__main__':
    find_crystal_classes()
