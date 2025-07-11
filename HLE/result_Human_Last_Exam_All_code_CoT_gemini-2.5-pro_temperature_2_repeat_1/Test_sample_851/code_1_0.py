def find_crystal_classes():
    """
    This function identifies crystal classes that are achiral, non-polar,
    and have the symmetry for optical activity (i.e., are non-centrosymmetric).
    """

    # Data for the 32 crystallographic point groups (crystal classes).
    # Each class is defined by its properties:
    # - chiral: True if it lacks any improper rotation axis (m or -n).
    # - polar: True if it belongs to one of the 10 polar classes.
    # - centrosymmetric: True if it possesses a center of inversion (-1).
    crystal_classes_data = {
        '1':    {'chiral': True,  'polar': True,  'centrosymmetric': False},
        '-1':   {'chiral': False, 'polar': False, 'centrosymmetric': True},
        '2':    {'chiral': True,  'polar': True,  'centrosymmetric': False},
        'm':    {'chiral': False, 'polar': True,  'centrosymmetric': False},
        '2/m':  {'chiral': False, 'polar': False, 'centrosymmetric': True},
        '222':  {'chiral': True,  'polar': False, 'centrosymmetric': False},
        'mm2':  {'chiral': False, 'polar': True,  'centrosymmetric': False},
        'mmm':  {'chiral': False, 'polar': False, 'centrosymmetric': True},
        '4':    {'chiral': True,  'polar': True,  'centrosymmetric': False},
        '-4':   {'chiral': False, 'polar': False, 'centrosymmetric': False},
        '4/m':  {'chiral': False, 'polar': False, 'centrosymmetric': True},
        '422':  {'chiral': True,  'polar': False, 'centrosymmetric': False},
        '4mm':  {'chiral': False, 'polar': True,  'centrosymmetric': False},
        '-42m': {'chiral': False, 'polar': False, 'centrosymmetric': False},
        '4/mmm':{'chiral': False, 'polar': False, 'centrosymmetric': True},
        '3':    {'chiral': True,  'polar': True,  'centrosymmetric': False},
        '-3':   {'chiral': False, 'polar': False, 'centrosymmetric': True},
        '32':   {'chiral': True,  'polar': False, 'centrosymmetric': False},
        '3m':   {'chiral': False, 'polar': True,  'centrosymmetric': False},
        '-3m':  {'chiral': False, 'polar': False, 'centrosymmetric': True},
        '6':    {'chiral': True,  'polar': True,  'centrosymmetric': False},
        '-6':   {'chiral': False, 'polar': False, 'centrosymmetric': False},
        '6/m':  {'chiral': False, 'polar': False, 'centrosymmetric': True},
        '622':  {'chiral': True,  'polar': False, 'centrosymmetric': False},
        '6mm':  {'chiral': False, 'polar': True,  'centrosymmetric': False},
        '-62m': {'chiral': False, 'polar': False, 'centrosymmetric': False},
        '6/mmm':{'chiral': False, 'polar': False, 'centrosymmetric': True},
        '23':   {'chiral': True,  'polar': False, 'centrosymmetric': False},
        'm-3':  {'chiral': False, 'polar': False, 'centrosymmetric': True},
        '432':  {'chiral': True,  'polar': False, 'centrosymmetric': False},
        '-43m': {'chiral': False, 'polar': False, 'centrosymmetric': False},
        'm-3m': {'chiral': False, 'polar': False, 'centrosymmetric': True}
    }

    print("Searching for crystal classes that meet the following three criteria:")
    print("1. Achiral (not chiral)")
    print("2. Non-polar (not polar)")
    print("3. Optically Active (not centrosymmetric)\n")

    result_classes = []
    for name, properties in crystal_classes_data.items():
        is_achiral = not properties['chiral']
        is_non_polar = not properties['polar']
        is_optically_active = not properties['centrosymmetric']

        if is_achiral and is_non_polar and is_optically_active:
            result_classes.append(name)
            
    # Sort for consistent output
    # Custom sort key to handle negative signs
    result_classes.sort(key=lambda x: (x.lstrip('-'), x.startswith('-')))

    print("The crystal classes satisfying all three conditions are:")
    # Using a loop to satisfy the instruction "output each number in the final equation!"
    for cls in result_classes:
      print(cls)

find_crystal_classes()