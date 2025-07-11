def find_crystal_classes():
    """
    Identifies achiral, non-polar crystal classes that are optically active.
    Hermann-Mauguin notation is used for the point groups.
    """
    # The 32 crystallographic point groups
    all_classes = {
        '1', '-1', '2', 'm', '2/m', '222', 'mm2', 'mmm', '4', '-4', '4/m',
        '422', '4mm', '-42m', '4/mmm', '3', '-3', '32', '3m', '-3m', '6',
        '-6', '6/m', '622', '6mm', '-6m2', '6/mmm', '23', 'm-3', '432',
        '-43m', 'm-3m'
    }

    # Condition 1: Optically active -> Must be non-centrosymmetric
    # The 11 centrosymmetric classes (possess a center of inversion) are not optically active.
    centrosymmetric_classes = {
        '-1', '2/m', 'mmm', '4/m', '4/mmm', '-3', '-3m', '6/m', '6/mmm', 'm-3', 'm-3m'
    }
    optically_active_classes = all_classes - centrosymmetric_classes

    # Condition 2: Achiral -> Must possess an improper rotation axis (m, -1, -3, -4, -6)
    # A class is achiral if it's not in the set of 11 purely chiral classes.
    chiral_classes = {
        '1', '2', '3', '4', '6', '222', '32', '422', '622', '23', '432'
    }
    achiral_classes = all_classes - chiral_classes

    # Condition 3: Non-polar -> Must not be in the set of 10 pyroelectric classes.
    polar_classes = {
        '1', '2', 'm', 'mm2', '3', '3m', '4', '4mm', '6', '6mm'
    }
    non_polar_classes = all_classes - polar_classes

    # Find the intersection of all three conditions
    result_set = optically_active_classes.intersection(achiral_classes).intersection(non_polar_classes)

    # Sort for a clean, ordered output
    result_list = sorted(list(result_set))

    print("The achiral and non-polar crystal classes that have the correct symmetry for optical activity are:")
    for cls in result_list:
        print(cls)

if __name__ == '__main__':
    find_crystal_classes()