def find_optically_active_achiral_nonpolar_classes():
    """
    Identifies crystal classes that are achiral, non-polar, and can exhibit optical activity.

    This is based on the following crystallographic principles:
    1.  Achiral: The class has an inversion center, a mirror plane, or a rotoinversion axis.
        These are any classes that are NOT in the chiral set.
    2.  Non-polar: The class is not one of the 10 polar classes.
    3.  Optically Active (Gyrotropic): The gyration tensor is non-zero. This requires the class
        to be non-centrosymmetric, but excludes a few special cases (3m, 4mm, 6mm, -62m) where
        the tensor is forced to be zero by other symmetries.
    """

    # Define properties for all 32 point groups
    chiral_classes = {'1', '2', '3', '4', '6', '222', '32', '422', '622', '23', '432'}
    polar_classes = {'1', '2', '3', '4', '6', 'm', 'mm2', '3m', '4mm', '6mm'}
    centrosymmetric_classes = {'-1', '2/m', 'mmm', '4/m', '-3', '6/m', '4/mmm', '-3m', '6/mmm', 'm-3', 'm-3m'}
    
    # Non-centrosymmetric classes where gyration tensor is identically zero
    zero_gyration_classes = {'3m', '4mm', '6mm', '-62m'}

    all_classes = {
        '1', '-1', '2', 'm', '2/m', '222', 'mm2', 'mmm', '4', '-4', '4/m',
        '422', '4mm', '-42m', '4/mmm', '3', '-3', '32', '3m', '-3m', '6',
        '-6', '6/m', '622', '6mm', '-62m', '6/mmm', '23', 'm-3', '432',
        '-43m', 'm-3m'
    }

    result_classes = []
    for cls in sorted(list(all_classes)):
        # Condition 1: Achiral (not in the chiral set)
        is_achiral = cls not in chiral_classes

        # Condition 2: Non-polar (not in the polar set)
        is_non_polar = cls not in polar_classes

        # Condition 3: Has correct symmetry for optical activity
        # This means it must be non-centrosymmetric AND not in the zero-gyration set.
        is_non_centrosymmetric = cls not in centrosymmetric_classes
        has_gyration = is_non_centrosymmetric and (cls not in zero_gyration_classes)

        if is_achiral and is_non_polar and has_gyration:
            result_classes.append(cls)

    print("The achiral, non-polar crystal classes that can exhibit optical activity are:")
    for cls in result_classes:
        print(cls)

find_optically_active_achiral_nonpolar_classes()