def find_special_crystal_classes():
    """
    Identifies crystal classes that are achiral, non-polar, and can exhibit optical activity.

    This is a fascinating case in crystallography where the common rule of thumb
    "optical activity requires chirality" is shown to have exceptions.

    The logic is as follows:
    1. A crystal class can exhibit optical activity if and only if it is non-centrosymmetric,
       with the special exception of classes -6 and -6m2, which are non-centrosymmetric but
       are still optically inactive due to other symmetry constraints.
    2. A crystal class is achiral if it is NOT chiral. A class is chiral if it lacks all
       improper rotation axes (inversion, mirror planes, rotoinversion).
    3. A crystal class is non-polar if it does not have a unique polar axis.

    This script defines the properties of all 32 point groups and filters them based on these
    three conditions.
    """

    # Data for all 32 point groups.
    # Format: Hermann-Mauguin Symbol: [is_centrosymmetric, is_chiral, is_polar]
    point_groups_data = {
        # Triclinic
        "1":     [False, True,  True],
        "-1":    [True,  False, False],
        # Monoclinic
        "2":     [False, True,  True],
        "m":     [False, False, True],
        "2/m":   [True,  False, False],
        # Orthorhombic
        "222":   [False, True,  False],
        "mm2":   [False, False, True],
        "mmm":   [True,  False, False],
        # Tetragonal
        "4":     [False, True,  True],
        "-4":    [False, False, False],
        "4/m":   [True,  False, False],
        "422":   [False, True,  False],
        "4mm":   [False, False, True],
        "-42m":  [False, False, False],
        "4/mmm": [True,  False, False],
        # Trigonal
        "3":     [False, True,  True],
        "-3":    [True,  False, False],
        "32":    [False, True,  False],
        "3m":    [False, False, True],
        "-3m":   [True,  False, False],
        # Hexagonal
        "6":     [False, True,  True],
        "-6":    [False, False, False],
        "6/m":   [True,  False, False],
        "622":   [False, True,  False],
        "6mm":   [False, False, True],
        "-6m2":  [False, False, False],
        "6/mmm": [True,  False, False],
        # Cubic
        "23":    [False, True,  False],
        "m-3":   [True,  False, False],
        "432":   [False, True,  False],
        "-43m":  [False, False, False],
        "m-3m":  [True,  False, False],
    }

    # These two non-centrosymmetric classes are known to be optically inactive.
    optically_inactive_exceptions = {"-6", "-6m2"}

    result_classes = []

    for name, properties in point_groups_data.items():
        is_centrosymmetric, is_chiral, is_polar = properties

        # Condition 1: Must be ACHIRAL (i.e., not chiral)
        if is_chiral:
            continue

        # Condition 2: Must be NON-POLAR
        if is_polar:
            continue

        # Condition 3: Must be OPTICALLY ACTIVE
        # This means it must be non-centrosymmetric AND not one of the special exceptions.
        is_optically_active = (not is_centrosymmetric) and (name not in optically_inactive_exceptions)

        if is_optically_active:
            result_classes.append(name)

    print("The achiral, non-polar crystal classes that can exhibit optical activity are:")
    # We join with a comma and space for readability, as requested
    print(", ".join(result_classes))

# Execute the function to find and print the answer.
find_special_crystal_classes()