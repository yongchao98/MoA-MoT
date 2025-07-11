def find_optically_active_classes():
    """
    Analyzes the 32 crystal classes to determine which ones fit certain criteria
    related to optical activity, chirality, and polarity.
    """

    # The 11 chiral (and therefore optically active) crystal classes.
    # They lack mirror planes and an inversion center.
    chiral_classes = {'1', '2', '222', '3', '32', '4', '422', '6', '622', '23', '432'}

    # The 21 achiral crystal classes (all other classes).
    all_32_classes = {
        '1', '-1', '2', 'm', '2/m', '222', 'mm2', 'mmm', '3', '-3', '32',
        '3m', '-3m', '4', '-4', '4/m', '422', '4mm', '-42m', '4/mmm', '6',
        '-6', '6/m', '622', '6mm', '-6m2', '6/mmm', '23', 'm-3', '432', '-43m', 'm-3m'
    }
    achiral_classes = all_32_classes - chiral_classes

    # The 10 polar (pyroelectric) crystal classes.
    # They possess a unique polar direction.
    polar_classes = {'1', '2', '3', '4', '6', 'm', 'mm2', '3m', '4mm', '6mm'}

    print("--- Analyzing the User's Original Question ---")
    print("Question: Which 'achiral' crystal classes are optically active?")
    print("Principle: Optical activity requires a crystal class to be 'chiral'.")
    print("Therefore, we are looking for classes that are in both the 'achiral' set and the 'chiral' set.")

    # Find the intersection of achiral and chiral sets. This is logically impossible.
    achiral_and_optically_active = achiral_classes.intersection(chiral_classes)

    print(f"\nResult: There are {len(achiral_and_optically_active)} such classes.")
    if not achiral_and_optically_active:
        print("This is because the condition is a contradiction, as explained.")

    print("\n--- Answering a Re-interpreted, More Likely Question ---")
    print("Re-interpreted Question: Which 'chiral' and 'non-polar' crystal classes are optically active?")
    print("Logic: We find the classes that are in the 'chiral' set but not in the 'polar' set.")

    # Optically active (chiral) and non-polar classes.
    chiral_non_polar_classes = chiral_classes.difference(polar_classes)

    print(f"\nResult: There are {len(chiral_non_polar_classes)} chiral and non-polar crystal classes.")
    print("These classes are optically active but not pyroelectric.")
    print("The classes are:")
    # Print each class name, which is a "number" in Hermann-Mauguin notation
    for cls in sorted(list(chiral_non_polar_classes)):
        print(cls)

if __name__ == '__main__':
    find_optically_active_classes()