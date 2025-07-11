def find_special_crystal_classes():
    """
    Analyzes the 32 crystallographic point groups to find those that are
    achiral and non-polar, yet can exhibit optical activity.

    A class can exhibit optical activity if it lacks a center of inversion and mirror planes.
    This script filters the 32 point groups based on these criteria.
    """
    # Data for the 32 point groups (Hermann-Mauguin notation)
    # Properties: is_chiral, is_polar, has_inversion_center, has_mirror_plane
    point_groups = [
        # Triclinic
        {'name': '1', 'is_chiral': True, 'is_polar': True, 'has_inversion_center': False, 'has_mirror_plane': False},
        {'name': '-1', 'is_chiral': False, 'is_polar': False, 'has_inversion_center': True, 'has_mirror_plane': False},
        # Monoclinic
        {'name': '2', 'is_chiral': True, 'is_polar': True, 'has_inversion_center': False, 'has_mirror_plane': False},
        {'name': 'm', 'is_chiral': False, 'is_polar': True, 'has_inversion_center': False, 'has_mirror_plane': True},
        {'name': '2/m', 'is_chiral': False, 'is_polar': False, 'has_inversion_center': True, 'has_mirror_plane': True},
        # Orthorhombic
        {'name': '222', 'is_chiral': True, 'is_polar': False, 'has_inversion_center': False, 'has_mirror_plane': False},
        {'name': 'mm2', 'is_chiral': False, 'is_polar': True, 'has_inversion_center': False, 'has_mirror_plane': True},
        {'name': 'mmm', 'is_chiral': False, 'is_polar': False, 'has_inversion_center': True, 'has_mirror_plane': True},
        # Tetragonal
        {'name': '4', 'is_chiral': True, 'is_polar': True, 'has_inversion_center': False, 'has_mirror_plane': False},
        {'name': '-4', 'is_chiral': False, 'is_polar': False, 'has_inversion_center': False, 'has_mirror_plane': False},
        {'name': '4/m', 'is_chiral': False, 'is_polar': False, 'has_inversion_center': True, 'has_mirror_plane': True},
        {'name': '422', 'is_chiral': True, 'is_polar': False, 'has_inversion_center': False, 'has_mirror_plane': False},
        {'name': '4mm', 'is_chiral': False, 'is_polar': True, 'has_inversion_center': False, 'has_mirror_plane': True},
        {'name': '-42m', 'is_chiral': False, 'is_polar': False, 'has_inversion_center': False, 'has_mirror_plane': True},
        {'name': '4/mmm', 'is_chiral': False, 'is_polar': False, 'has_inversion_center': True, 'has_mirror_plane': True},
        # Trigonal
        {'name': '3', 'is_chiral': True, 'is_polar': True, 'has_inversion_center': False, 'has_mirror_plane': False},
        {'name': '-3', 'is_chiral': False, 'is_polar': False, 'has_inversion_center': True, 'has_mirror_plane': False},
        {'name': '32', 'is_chiral': True, 'is_polar': False, 'has_inversion_center': False, 'has_mirror_plane': False},
        {'name': '3m', 'is_chiral': False, 'is_polar': True, 'has_inversion_center': False, 'has_mirror_plane': True},
        {'name': '-3m', 'is_chiral': False, 'is_polar': False, 'has_inversion_center': True, 'has_mirror_plane': True},
        # Hexagonal
        {'name': '6', 'is_chiral': True, 'is_polar': True, 'has_inversion_center': False, 'has_mirror_plane': False},
        {'name': '-6', 'is_chiral': False, 'is_polar': False, 'has_inversion_center': False, 'has_mirror_plane': True},
        {'name': '6/m', 'is_chiral': False, 'is_polar': False, 'has_inversion_center': True, 'has_mirror_plane': True},
        {'name': '622', 'is_chiral': True, 'is_polar': False, 'has_inversion_center': False, 'has_mirror_plane': False},
        {'name': '6mm', 'is_chiral': False, 'is_polar': True, 'has_inversion_center': False, 'has_mirror_plane': True},
        {'name': '-6m2', 'is_chiral': False, 'is_polar': False, 'has_inversion_center': False, 'has_mirror_plane': True},
        {'name': '6/mmm', 'is_chiral': False, 'is_polar': False, 'has_inversion_center': True, 'has_mirror_plane': True},
        # Cubic
        {'name': '23', 'is_chiral': True, 'is_polar': False, 'has_inversion_center': False, 'has_mirror_plane': False},
        {'name': 'm-3', 'is_chiral': False, 'is_polar': False, 'has_inversion_center': True, 'has_mirror_plane': True},
        {'name': '432', 'is_chiral': True, 'is_polar': False, 'has_inversion_center': False, 'has_mirror_plane': False},
        {'name': '-43m', 'is_chiral': False, 'is_polar': False, 'has_inversion_center': False, 'has_mirror_plane': True},
        {'name': 'm-3m', 'is_chiral': False, 'is_polar': False, 'has_inversion_center': True, 'has_mirror_plane': True},
    ]

    print("Searching for crystal classes with the following properties:")
    print("1. Achiral (is_chiral = False)")
    print("2. Non-polar (is_polar = False)")
    print("3. Can be optically active (has_inversion_center = False AND has_mirror_plane = False)")
    print("-" * 30)

    results = []
    for pg in point_groups:
        # Condition 1: Achiral
        is_achiral = not pg['is_chiral']
        # Condition 2: Non-polar
        is_non_polar = not pg['is_polar']
        # Condition 3: Can be optically active (no inversion, no mirror)
        can_be_active = not pg['has_inversion_center'] and not pg['has_mirror_plane']

        if is_achiral and is_non_polar and can_be_active:
            results.append(pg['name'])

    if results:
        print("Found matching crystal class(es):")
        for r in results:
            # We must print the Hermann-Mauguin symbol character by character
            # For '-4' it would be print('-', '4')
            if r.startswith('-'):
                 print("The crystal class is:", r[0], r[1:])
            else:
                 print("The crystal class is:", r)
    else:
        print("No crystal classes found matching all criteria.")

find_special_crystal_classes()