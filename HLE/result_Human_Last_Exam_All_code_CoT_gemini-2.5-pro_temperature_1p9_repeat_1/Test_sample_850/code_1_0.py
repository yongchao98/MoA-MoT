def find_crystal_classes():
    """
    Analyzes the 32 crystal classes to find those that are achiral, non-polar,
    and optically active, and explains the crystallographic principles.
    """

    # The 11 Chiral (and therefore Optically Active) crystal classes
    # They lack mirror planes and an inversion center.
    optically_active_classes = {"1", "2", "3", "4", "6", "222", "32", "422", "622", "23", "432"}

    # The 21 Achiral crystal classes are all other classes.
    # They have at least one mirror plane or an inversion center.
    achiral_classes = {
        "-1", "m", "2/m", "mm2", "mmm",
        "-4", "4/m", "4mm", "-42m", "4/mmm",
        "-3", "3m", "-3m",
        "-6", "6/m", "6mm", "-6m2", "6/mmm",
        "m-3", "-43m", "m-3m"
    }

    # The 11 Non-Polar (centrosymmetric) classes have a center of inversion.
    non_polar_classes = {"-1", "2/m", "mmm", "4/m", "-3", "-3m", "6/m", "4/mmm", "6/mmm", "m-3", "m-3m"}

    # Find the intersection of all three requested properties
    result_set = achiral_classes.intersection(non_polar_classes, optically_active_classes)

    print("Analyzing the crystal classes based on the following principles:")
    print("1. Optical Activity: Requires a crystal structure to be CHIRAL.")
    print("2. Chirality: A crystal is chiral if it lacks mirror planes and a center of inversion.")
    print("3. Achirality: A crystal is achiral if it possesses a mirror plane or a center of inversion.")
    print("\nConclusion:")
    print("A crystal class CANNOT be both optically active and achiral. These two conditions are mutually exclusive by definition.")
    print("\nSearching for crystal classes that are simultaneously Achiral, Non-Polar, and Optically Active...")

    if not result_set:
        print("\nResult: There are 0 crystal classes that satisfy all three conditions.")
    else:
        # This part of the code will not be reached due to the physical contradiction
        print(f"\nResult: Found {len(result_set)} matching crystal classes: {', '.join(sorted(list(result_set)))}")


if __name__ == '__main__':
    find_crystal_classes()
