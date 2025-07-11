import textwrap

def find_contradictory_crystal_classes():
    """
    This script identifies crystal classes based on symmetry properties to find
    any that are achiral, non-polar, and optically active simultaneously.
    """

    # According to crystallography, a crystal class is optically active if and only if it is chiral.
    # These classes lack a center of inversion and mirror planes.
    optically_active_classes = {
        '1', '2', '3', '4', '6', '222', '32', '422', '622', '23', '432'
    }

    # A class is non-polar if it is centrosymmetric (possesses a center of inversion).
    non_polar_classes = {
        '-1', '2/m', 'mmm', '-3', '-3m', '4/m', '4/mmm', '6/m', '6/mmm', 'm-3', 'm-3m'
    }

    # A class is achiral if it is NOT chiral (i.e., it possesses a center of inversion or a mirror plane).
    # We can define all 32 classes and subtract the chiral (optically active) ones.
    all_32_crystal_classes = {
        '1', '-1', '2', 'm', '2/m', '222', 'mm2', 'mmm',
        '4', '-4', '4/m', '422', '4mm', '-42m', '4/mmm',
        '3', '-3', '32', '3m', '-3m',
        '6', '-6', '6/m', '622', '6mm', '-6m2', '6/mmm',
        '23', 'm-3', '432', '-43m', 'm-3m'
    }
    achiral_classes = all_32_crystal_classes - optically_active_classes

    # --- Analysis ---
    print("Step 1: Define the required sets of crystal classes based on symmetry.")
    print("-" * 70)

    print("Set A: Achiral Classes (must have mirror or inversion center)")
    print(textwrap.fill(str(sorted(list(achiral_classes))), width=70))
    print("\nSet B: Non-polar Classes (must have inversion center)")
    print(textwrap.fill(str(sorted(list(non_polar_classes))), width=70))
    print("\nSet C: Optically Active Classes (must NOT have mirror or inversion center)")
    print(textwrap.fill(str(sorted(list(optically_active_classes))), width=70))
    print("-" * 70)

    # Find the intersection of all three sets.
    # The equation is: Result = A \u2229 B \u2229 C
    result = achiral_classes.intersection(non_polar_classes).intersection(optically_active_classes)

    print("\nStep 2: Find the intersection of all three sets.")
    print("Equation: Result = (Achiral \u2229 Non-polar \u2229 Optically Active)")
    print("-" * 70)
    print(f"Number of classes in Achiral set: {len(achiral_classes)}")
    print(f"Number of classes in Non-polar set: {len(non_polar_classes)}")
    print(f"Number of classes in Optically Active set: {len(optically_active_classes)}")
    print("\nFinal Result of Intersection:")
    if not result:
        print("The resulting set is empty: {}")
    else:
        print(result)

    print("-" * 70)
    print("\nConclusion:")
    conclusion_text = (
        "There are no crystal classes that are simultaneously achiral, non-polar, "
        "and optically active. The conditions are mutually exclusive. A class must "
        "lack an inversion center to be optically active, but it must possess an "
        "inversion center to be non-polar."
    )
    print(textwrap.fill(conclusion_text, width=70))

if __name__ == "__main__":
    find_contradictory_crystal_classes()