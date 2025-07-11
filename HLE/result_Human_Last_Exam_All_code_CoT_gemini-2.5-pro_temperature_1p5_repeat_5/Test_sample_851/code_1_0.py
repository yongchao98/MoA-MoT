def find_crystal_classes():
    """
    This script identifies crystal classes that are achiral, non-polar, and optically active.
    """
    # The 32 crystallographic point groups (Hermann-Mauguin notation)
    all_classes = {'1', '-1', '2', 'm', '2/m', '222', 'mm2', 'mmm', '4', '-4', '4/m', '422', '4mm', '-42m', '4/mmm', '3', '-3', '32', '3m', '-3m', '6', '-6', '6/m', '622', '6mm', '-62m', '6/mmm', '23', 'm-3', '432', '-43m', 'm-3m'}

    # 1. CRITERION: OPTICAL ACTIVITY
    # Optically inactive classes include all 11 centrosymmetric classes plus the special case -43m.
    optically_inactive = {'-1', '2/m', 'mmm', '4/m', '4/mmm', '-3', '-3m', '6/m', '6/mmm', 'm-3', 'm-3m', '-43m'}
    optically_active = all_classes - optically_inactive

    # 2. CRITERION: ACHIRALITY
    # Chiral classes are those that do not have any improper rotation axes (-n, including mirror planes).
    # Achiral classes are all other classes.
    chiral = {'1', '2', '222', '3', '32', '4', '422', '6', '622', '23', '432'}
    achiral = all_classes - chiral

    # 3. CRITERION: NON-POLARITY
    # Polar classes are those with a unique polar direction (pyroelectric classes).
    # Non-polar classes are all other classes.
    polar = {'1', '2', '3', '4', '6', 'm', 'mm2', '3m', '4mm', '6mm'}
    non_polar = all_classes - polar

    # Find the intersection of all three criteria
    result_classes = optically_active.intersection(achiral).intersection(non_polar)

    # Sort the results for consistent output
    final_list = sorted(list(result_classes))

    print("Step 1: The crystal classes must be optically active.")
    print(f"There are {len(optically_active)} such classes: {sorted(list(optically_active))}\n")

    print("Step 2: From this group, we select the ones that are achiral.")
    achiral_and_active = optically_active.intersection(achiral)
    print(f"This leaves {len(achiral_and_active)} classes: {sorted(list(achiral_and_active))}\n")
    
    print("Step 3: From this smaller group, we select the ones that are non-polar.")
    print(f"This leaves the final set of {len(final_list)} classes.\n")

    print("The final list of crystal classes that are achiral, non-polar, and optically active is:")
    # The prompt requests to "output each number in the final equation".
    # This is interpreted as printing the final list of classes.
    for crystal_class in final_list:
        print(crystal_class)

find_crystal_classes()