import textwrap

def analyze_crystal_classes():
    """
    Analyzes the 32 crystal classes to find those that are achiral,
    non-polar, and optically active, and explains the result.
    """
    # Define all 32 crystallographic point groups (crystal classes)
    all_classes = {
        '1', 'ī', '2', 'm', '2/m', '222', 'mm2', 'mmm', '4', '4̄', '4/m',
        '422', '4mm', '4̄2m', '4/mmm', '3', '3̄', '32', '3m', '3̄m', '6',
        '6̄', '6/m', '622', '6mm', '6̄m2', '6/mmm', '23', 'm3̄', '432',
        '4̄3m', 'm3̄m'
    }

    # Define subsets based on physical properties and symmetry.

    # 1. Optical Activity:
    # A crystal is optically active if and only if it is chiral.
    # Chiral classes lack any roto-inversion axes (i.e., no mirror planes 'm'
    # or a center of inversion 'ī').
    optically_active_classes = {
        '1', '2', '3', '4', '6', '222', '32', '422', '622', '23', '432'
    }

    # 2. Achirality:
    # Achiral classes are all classes that are NOT chiral.
    achiral_classes = all_classes - optically_active_classes

    # 3. Polarity:
    # Polar classes possess a unique vector direction (a polar axis).
    polar_classes = {'1', '2', 'm', 'mm2', '3', '3m', '4', '4mm', '6', '6mm'}
    # Non-polar classes are all classes that are NOT polar.
    non_polar_classes = all_classes - polar_classes

    # --- Explanation ---
    print("Analysis of Crystal Classes by Symmetry Properties")
    print("=" * 70)
    explanation = (
        "The query asks for crystal classes that are simultaneously achiral, non-polar, "
        "and have the correct symmetry for optical activity. A fundamental principle "
        "of physics is that intrinsic optical activity is only possible in chiral structures. "
        "A chiral structure is one that is not superimposable on its mirror image. "
        "Therefore, by definition, a crystal must be CHIRAL to be optically active. "
        "The request for an ACHIRAL optically active crystal class is a contradiction."
    )
    print("\n".join(textwrap.wrap(explanation, 70)))
    print("\nTo demonstrate this, we will find the intersection of the following sets:")

    # --- Set Definitions ---
    print("\n1. Optically Active (i.e., Chiral) Classes:")
    print(f"   - Total: {len(optically_active_classes)}")
    print(f"   - Members: {', '.join(sorted(list(optically_active_classes)))}")

    print("\n2. Achiral Classes:")
    print(f"   - Total: {len(achiral_classes)}")
    print(f"   - Members: {', '.join(sorted(list(achiral_classes)))}")

    print("\n3. Non-Polar Classes:")
    print(f"   - Total: {len(non_polar_classes)}")
    print(f"   - Members: {', '.join(sorted(list(non_polar_classes)))}")

    # --- Calculation ---
    print("\n" + "=" * 70)
    print("Calculating the intersection of {Achiral} AND {Non-Polar} AND {Optically Active} sets...")

    # Find the intersection of the three sets
    result_set = achiral_classes.intersection(non_polar_classes).intersection(optically_active_classes)

    # --- Result ---
    print("\n--- FINAL RESULT ---")
    if not result_set:
        print("The resulting set is empty.")
        print("\nConclusion: There are NO crystal classes that satisfy all three conditions.")
    else:
        # This code block should not be reachable due to the physical principles.
        print(f"Found {len(result_set)} matching crystal classes:")
        print(', '.join(sorted(list(result_set))))

if __name__ == '__main__':
    analyze_crystal_classes()